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
  //case class PopZoom(n:Int) extends TracePlotterMsg
  //case class SetZoomFilter(f:(Datum, Int) => Boolean) extends TracePlotterMsg
  case object HideLegend extends TracePlotterMsg
  case object ShowLegend extends TracePlotterMsg
  case class RemoveAssayTracePeak(peakID: Int) extends TracePlotterMsg

  case class PlotUpdate(id: PlotID, plot: BufferedImage, control: PlotsControl[Datum, Datum, Datum])
}

class TracePlotter(swingActor: ActorRef, id: PlotID, hideLegend: Boolean) extends Actor {

  import TracePlotter._
  import Scale._
  import Stratifier._
  import PlotBuffer._
  import PeakIntegrator._

  var size = new Dimension(1000, 600)
  var plot: Option[LinePlot[Datum]] = None
  val assayPeakPlot = new AssayPeakPlot
  val peakIDs = new IDGenerator

  def receive = {
    case (reply: MasterReply, assay: Assay) =>
      if (reply.hasStatus)
        println(reply.getStatus.getStatusMsg)
      if (reply.hasTraces) {
        plot = Some(getTracePlot(reply.getTraces, assay))
        swingActor ! plotUpdate(plot.get)
      }

    case msg: MSDataProtocolActors.MSDataProtocolMsg =>
      println(msg)

    case NewSize(_, d) =>
      if (d.getHeight > 0 && d.getWidth > 0)
        size = d
      plot.foreach(p => swingActor ! plotUpdate(p))

    case Select(fromID, f) =>
      println("NEW SELECT from %s to %s ".format(fromID, id))
      for (p <- plot) {
        val integrator = context.actorOf(Props[PeakIntegrator], "integrator")
        def filter(data: Seq[Datum], f: (Datum, Int) => Boolean) =
          for {
            i <- 0 until data.length
            if f(data(i), i)
          } yield data(i)

        integrator ! PeakIntegrator.Integrate(id, peakIDs.next, filter(p.data, f), swingActor)
      }

    case atp: AssayTracePeak =>
      assayPeakPlot.peaks.clear
      assayPeakPlot.peaks += atp
      for (p <- plot)
        swingActor ! plotUpdate(p)

    case RemoveAssayTracePeak(peakID) =>
      assayPeakPlot.remove(peakID)
      for (p <- plot)
        swingActor ! plotUpdate(p)

    case PopZoom(_, n) =>
      for (p <- plot)
        if (p.filters.nonEmpty) {
          p.filters.pop
          swingActor ! plotUpdate(p)
        }

    case NewZoom(_, f) =>
      for (p <- plot) {
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

  def getTracePlot(t: Traces, assay: Assay) = {
    val data = new ArrayBuffer[Datum]
    for (prec <- assay.precs) {
      t.getPrecursorList.find(t =>
        prec.mz > t.getPrecursor.getLmz &&
          prec.mz < t.getPrecursor.getHmz) match {
        case None => println("TracePlotter: cannot match trace with assay!")
        case Some(pTrace) =>
          val trace = pTrace.getTrace.getTimeList.zip(pTrace.getTrace.getIntensityList)
          data ++= trace.map(t => Datum(t._1, t._2, prec.id))
      }
    }

    for (i <- 0 until assay.frags.size) {
      val temp = assay.frags(i)
      t.getFragmentList.find(t =>
        temp.precMz > t.getFragment.getPrecursor.getLmz &&
          temp.precMz < t.getFragment.getPrecursor.getHmz &&
          temp.fragMz > t.getFragment.getFragment.getLmz &&
          temp.fragMz < t.getFragment.getFragment.getHmz) match {
        case None => println("TracePlotter: cannot match trace with assay!")
        case Some(fTrace) =>
          val temp1 = fTrace.getTrace.getIntensityList
          val temp2 = savitzkyGolay9(temp1.asScala.map(_.doubleValue).toSeq)
          val temp3 = testWavelet(temp1.asScala.map(_.doubleValue).toSeq, 0)
          val temp4 = testWavelet(temp1.asScala.map(_.doubleValue).toSeq, 1)

          val signalDB3 = daubchies3(temp1.asScala.map(_.doubleValue).toSeq, 0)
          val signalWL0 = testWavelet(temp1.asScala.map(_.doubleValue).toSeq, 0)
          val signalDB8 = daubchies8(temp1.asScala.map(_.doubleValue).toSeq, 0)
          val trace1 = fTrace.getTrace.getTimeList.zip(temp1)
//          val trace2 = fTrace.getTrace.getTimeList.zip(signalWL0)
//          val trace3 = fTrace.getTrace.getTimeList.zip(signalDB3)
//          val trace4 = fTrace.getTrace.getTimeList.zip(signalDB8)
          data ++= trace1.map(t => Datum(t._1, t._2, temp.id + " org"))
//          data ++= trace2.map(t => Datum(t._1, t._2, temp.id + " db2"))
//          data ++= trace3.map(t => Datum(t._1, t._2, temp.id + " db3"))
//
//          data ++= trace4.map(t => Datum(t._1, t._2, temp.id + " db8"))
      }
    }
    val p = new LinePlot(data)
      .x(_.rt)
      .y(_.intensity)
      .color(_.id)
    p.hideLegends = hideLegend
    p.overlays += assayPeakPlot
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
  def testWavelet(x: Seq[Double], levels: Int): Seq[Double] = {
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
    if (levels > 0) {
      yl = testWavelet(yl.toSeq, levels - 1).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 1) % yl.size + yl.size) % yl.size) * H(1) + yl((i) % yl.size) * H(3) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G(1) + yh((yh.size + i) % yh.size) * G(3)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H(0) + yl((yl.size + i) % yl.size) * H(2) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G(0) + yh((yh.size + i) % yh.size) * G(2)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
  def daubchies3(x: Seq[Double], levels: Int): Seq[Double] = {

    val H3 = Array(
      math.sqrt(2) * (1 + math.sqrt(10) + math.sqrt(5 + 2 * math.sqrt(10))) / 32,
      math.sqrt(2) * (5 + math.sqrt(10) + 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
      math.sqrt(2) * (10 - 2 * math.sqrt(10) + 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
      math.sqrt(2) * (10 - 2 * math.sqrt(10) - 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
      math.sqrt(2) * (5 + math.sqrt(10) - 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
      math.sqrt(2) * (1 + math.sqrt(10) - math.sqrt(5 + 2 * math.sqrt(10))) / 32)
    val G3 = Array(H3(5), -H3(4), H3(3), -H3(2), H3(1), -H3(0))
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until 2) {
      yh(i) = 0
      yl(i) = x(i)
    }
    for (i <- 2 until x.size / 2) //forward Transform
    {
      yh(i) = 0
      //    yh(i)= x((i*2)%x.size)*G2(3) + x((i*2+1)%x.size)*G2(2) + x((i*2+2)%x.size)*G2(1) + x((i*2+3)%x.size)*G2(0)
      yl(i) = x((i * 2) % x.size) * H3(5) + x((i * 2 + 1) % x.size) * H3(4) + x((i * 2 + 2) % x.size) * H3(3) + x((i * 2 + 3) % x.size) * H3(2) + x((i * 2 + 4) % x.size) * H3(1) + x((i * 2 + 5) % x.size) * H3(0)
    }
    if (levels > 0) {
      yl = daubchies3(yl.toSeq, levels - 1).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 2) % yl.size + yl.size) % yl.size) * H3(1) + yl(((i - 1) % yl.size + yl.size) % yl.size) * H3(3) + yl((i) % yl.size) * H3(5) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G3(1) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G3(3) + yh((yh.size + i) % yh.size) * G3(5)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 2 + i) % yl.size + yl.size) % yl.size) * H3(0) + yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H3(2) + yl((yl.size + i) % yl.size) * H3(4) + yh(((yh.size - 2 + i) % yh.size + yh.size) % yh.size) * G3(0) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G3(2) + yh((yh.size + i) % yh.size) * G3(4)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
  def daubchies8(x: Seq[Double], levels: Int): Seq[Double] = {
    val H8 = Array(
      0.0018899503329007,
      0.0003029205145516,
      -0.0149522583367926,
      -0.0038087520140601,
      0.0491371796734768,
      0.0272190299168137,
      -0.0519458381078751,
      -0.3644418948359564,
      0.7771857516997478,
      -0.4813596512592012,
      -0.0612733590679088,
      0.1432942383510542,
      0.0076074873252848,
      -0.0316950878103452,
      -0.0005421323316355,
      0.0033824159513594)
    val G8 = Array(
      -0.0033824159513594,
      -0.0005421323316355,
      0.0316950878103452,
      0.0076074873252848,
      -0.1432942383510542,
      -0.0612733590679088,
      0.4813596512592012,
      0.7771857516997478,
      0.3644418948359564,
      -0.0519458381078751,
      -0.0272190299168137,
      0.0491371796734768,
      0.0038087520140601,
      -0.0149522583367926,
      -0.0003029205145516,
      0.0018899503329007)
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform
    {
      
      
      yh(i) = x((i * 2) % x.size) * G8(15) + x((i * 2 + 1) % x.size) * G8(14) + x((i * 2 + 2) % x.size) * G8(13) + x((i * 2 + 3) % x.size) * G8(12) + x((i * 2 + 4) % x.size) * G8(11) + x((i * 2 + 5) % x.size) * G8(10) + x((i * 2 + 6) % x.size) * G8(9) + x((i * 2 + 7) % x.size) * G8(8) + x((i * 2 + 8) % x.size) * G8(7) + x((i * 2 + 9) % x.size) * G8(6) + x((i * 2 + 10) % x.size) * G8(5) + x((i * 2 + 11) % x.size) * G8(4) + x((i * 2 + 12) % x.size) * G8(3) + x((i * 2 + 13) % x.size) * G8(2) + x((i * 2 + 14) % x.size) * G8(1) + x((i * 2 + 15) % x.size) * G8(0)
      if(yh(i)<50){
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * H8(15) + x((i * 2 + 1) % x.size) * H8(14) + x((i * 2 + 2) % x.size) * H8(13) + x((i * 2 + 3) % x.size) * H8(12) + x((i * 2 + 4) % x.size) * H8(11) + x((i * 2 + 5) % x.size) * H8(10) + x((i * 2 + 6) % x.size) * H8(9) + x((i * 2 + 7) % x.size) * H8(8) + x((i * 2 + 8) % x.size) * H8(7) + x((i * 2 + 9) % x.size) * H8(6) + x((i * 2 + 10) % x.size) * H8(5) + x((i * 2 + 11) % x.size) * H8(4) + x((i * 2 + 12) % x.size) * H8(3) + x((i * 2 + 13) % x.size) * H8(2) + x((i * 2 + 14) % x.size) * H8(1) + x((i * 2 + 15) % x.size) * H8(0)
    }
    if (levels > 0) {
      yl = daubchies3(yl.toSeq, levels - 1).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 7) % yl.size + yl.size) % yl.size) * H8(1) + yl(((i - 6) % yl.size + yl.size) % yl.size) * H8(3) + yl(((i - 5) % yl.size + yl.size) % yl.size) * H8(5) + yl(((i - 4) % yl.size + yl.size) % yl.size) * H8(7) + yl(((i - 3) % yl.size + yl.size) % yl.size) * H8(9) + yl(((i - 2) % yl.size + yl.size) % yl.size) * H8(11) + yl(((i - 1) % yl.size + yl.size) % yl.size) * H8(13) + yl(((i) % yl.size + yl.size) % yl.size) * H8(15) + yh((((yh.size - 7 + i) % yh.size) + yh.size) % yh.size) * G8(1) + yh((((yh.size - 6 + i) % yh.size) + yh.size) % yh.size) * G8(3) + yh((((yh.size - 5 + i) % yh.size) + yh.size) % yh.size) * G8(5) + yh((((yh.size - 4 + i) % yh.size) + yh.size) % yh.size) * G8(7) + yh((((yh.size - 3 + i) % yh.size) + yh.size) % yh.size) * G8(9) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G8(11) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G8(13) + yh((((yh.size + i) % yh.size) + yh.size) % yh.size) * G8(15)
      temp((i * 2 + 1) % x.size) = yl(((yl.size + i - 7) % yl.size + yl.size) % yl.size) * H8(0) + yl(((yl.size + i - 6) % yl.size + yl.size) % yl.size) * H8(2) + yl(((yl.size + i - 5) % yl.size + yl.size) % yl.size) * H8(4) + yl(((yl.size + i - 4) % yl.size + yl.size) % yl.size) * H8(6) + yl(((yl.size + i - 3) % yl.size + yl.size) % yl.size) * H8(8) + yl(((yl.size + i - 2) % yl.size + yl.size) % yl.size) * H8(10) + yl(((yl.size + i - 1) % yl.size + yl.size) % yl.size) * H8(12) + yl(((yl.size + i) % yl.size + yl.size) % yl.size) * H8(14) + yh((((yh.size - 7 + i) % yh.size) + yh.size) % yh.size) * G8(1) + yh((((yh.size - 6 + i) % yh.size) + yh.size) % yh.size) * G8(3) + yh((((yh.size - 5 + i) % yh.size) + yh.size) % yh.size) * G8(5) + yh((((yh.size - 4 + i) % yh.size) + yh.size) % yh.size) * G8(7) + yh((((yh.size - 3 + i) % yh.size) + yh.size) % yh.size) * G8(9) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G8(11) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G8(13) + yh((((yh.size + i) % yh.size) + yh.size) % yh.size) * G8(15)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
}