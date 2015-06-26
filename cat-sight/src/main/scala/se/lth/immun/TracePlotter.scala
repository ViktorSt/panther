package se.lth.immun

import akka.actor._
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConversions._
import scala.swing.Dimension

import se.lth.immun.protocol.MSDataProtocolActors
import se.lth.immun.protocol.MSDataProtocol.Traces
import se.lth.immun.protocol.MSDataProtocol.MasterReply

import se.jt.{ Scale, Stratifier, LinePlot, Util, PlotsControl }
import java.awt.image.BufferedImage

import scala.collection.JavaConversions._

import collection.JavaConverters._
import CatSightPrimaries._

object TracePlotter {
  def props(swingActor: ActorRef, id: PlotID, hideLegend: Boolean) =
    Props(classOf[TracePlotter], swingActor, id, hideLegend)

  case class Datum(rt: Double, intensity: Double, id: String)
  trait TracePlotterMsg
  case class PopZoom(n: Int) extends TracePlotterMsg
  case class SetZoomFilter(f: (Datum, Int) => Boolean) extends TracePlotterMsg
  case object HideLegend extends TracePlotterMsg
  case object ShowLegend extends TracePlotterMsg

  case class PlotUpdate(id: PlotID, plot: BufferedImage, control: PlotsControl[Datum, Datum, Datum])
}

class TracePlotter(swingActor: ActorRef, id: PlotID, hideLegend: Boolean) extends Actor {

  import TracePlotter._
  import Scale._
  import Stratifier._

  var size = new Dimension(1000, 600)
  var plot: Option[LinePlot[Datum]] = None

  def receive = {
    case reply: MasterReply =>
      if (reply.hasStatus)
        println(reply.getStatus.getStatusMsg)
      if (reply.hasTraces) {
        plot = Some(getTracePlot(reply.getTraces))
        swingActor ! plotUpdate(plot.get)
      }

    case msg: MSDataProtocolActors.MSDataProtocolMsg =>
      println(msg)

    case d: Dimension =>
      if (d.getHeight > 0 && d.getWidth > 0)
        size = d
      plot.foreach(p => swingActor ! plotUpdate(p))

    case PopZoom(n) =>
      plot match {
        case Some(p) if p.filters.nonEmpty =>
          p.filters.pop
          swingActor ! plotUpdate(p)

        case Some(_) => {}
        case None    => {}
      }

    case SetZoomFilter(f) =>
      plot match {
        case None => {}
        case Some(p) =>
          p.filters.push(f)
          swingActor ! plotUpdate(p)
      }

    case HideLegend =>
      for (p <- plot)
        if (!p.hideLegends) {
          p.hideLegends = true
          swingActor ! plotUpdate(p)
        }

    case ShowLegend =>
      for (p <- plot)
        if (p.hideLegends) {
          p.hideLegends = false
          swingActor ! plotUpdate(p)
        }
  }

  def plotUpdate(p: LinePlot[Datum]) = {
    val imgCtrl = Util.drawToBuffer(size.width, size.height, p)
    PlotUpdate(id, imgCtrl.img, imgCtrl.control)
  }

  def getTracePlot(t: Traces) = {
    val data = new ArrayBuffer[Datum]
//    for (pTrace <- t.getPrecursorList) {
//      val id = "prec %.4f".format(pTrace.getPrecursor.getLmz)
//      val trace = pTrace.getTrace.getTimeList.zip(pTrace.getTrace.getIntensityList)
//      data ++= trace.map(t => Datum(t._1, t._2, id))
//    }
    for (i<-0 until t.getFragmentList.size) {
      val fTrace = t.getFragmentList()(i)
      val id = "frag %.4f -> %.4f".format(fTrace.getFragment.getPrecursor.getLmz, fTrace.getFragment.getFragment.getLmz)
      val temp1 = fTrace.getTrace.getIntensityList
      val temp2 = savitzkyGolay9(fTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
      val temp3 = testWavelet(fTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
      val trace1 = fTrace.getTrace.getTimeList.zip(temp1)
      val trace2 = fTrace.getTrace.getTimeList.zip(temp2)
      val trace3 = fTrace.getTrace.getTimeList.zip(temp3)
      data ++= trace1.map(t => Datum(t._1, t._2, id+ " org"))
      data ++= trace2.map(t => Datum(t._1, t._2, id+ " savitzky"))
      data ++= trace3.map(t => Datum(t._1, t._2, id+ " wavelet"))
    }
    val p = new LinePlot(data)
      .x(_.rt)
      .y(_.intensity)
      .color(_.id)
    p.hideLegends = hideLegend
    p
  }
  def savitzkyGolay9(a: Seq[Double]): Seq[Double] = {

    if (a.length < 9)
      return a

    var retArr = new Array[Double](a.length)
    var j = a.length
    for (i <- 0 until 4) {
      j -= 1
      retArr(i) = a(i)
      retArr(j) = a(j)
    }

    for (i <- 4 until a.length - 4) {
      retArr(i) = (0.417 * a(i) +
        +0.315 * (a(i - 1) + a(i + 1))
        + 0.070 * (a(i - 2) + a(i + 2))
        - 0.128 * (a(i - 3) + a(i + 3))
        + 0.035 * (a(i - 4) + a(i + 4)))
    }
    return retArr map (d => math.max(0, d))
  }
  def testWavelet(x: Seq[Double]): Seq[Double] = {
    val H = Array(
      (1 + math.sqrt(3)) / (4 * math.sqrt(2)),
      (3 + math.sqrt(3)) / (4 * math.sqrt(2)),
      (3 - math.sqrt(3)) / (4 * math.sqrt(2)),
      (1 - math.sqrt(3)) / (4 * math.sqrt(2)))

    val G = Array(H(3), -H(2), H(1), -H(0))
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform
    {
      yh(i) = 0
      //    yh(i)= x((i*2)%x.size)*G(3) + x((i*2+1)%x.size)*G(2) + x((i*2+2)%x.size)*G(1) + x((i*2+3)%x.size)*G(0)
      yl(i) = x((i * 2) % x.size) * H(3) + x((i * 2 + 1) % x.size) * H(2) + x((i * 2 + 2) % x.size) * H(1) + x((i * 2 + 3) % x.size) * H(0)
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 1) % yl.size + yl.size) % yl.size) * H(1) + yl((i) % yl.size) * H(3) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G(1) + yh((yh.size + i) % yh.size) * G(3)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H(0) + yl((yl.size + i) % yl.size) * H(2) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G(0) + yh((yh.size + i) % yh.size) * G(2)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
}