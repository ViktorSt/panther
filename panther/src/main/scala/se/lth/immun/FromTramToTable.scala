package se.lth.immun
import java.io.File
import java.io.FileReader
import java.io.BufferedReader
import akka.actor._
import akka.actor.{ Actor, ActorRef, Props }
import java.net.InetSocketAddress
import se.lth.immun.protocol.MSDataProtocolActors
import se.lth.immun.protocol.MSDataProtocol._
import se.lth.immun.xml.XmlReader
import java.util.Scanner
import scala.collection.mutable.ArrayBuffer
import scala.collection.{ JavaConversions => JConversions }
import scala.collection.mutable.Queue
import se.lth.immun.traml.ghost.GhostTraML
import scala.collection.JavaConversions._
import collection.JavaConverters._
import akka.io.Tcp
import Tcp._
import se.lth.immun.DataModel._
import se.lth.immun.Filters._
import se.lth.immun.DetectionAndGrouping._
import se.lth.immun.PeakFinder._
import se.lth.immun.traml.Peptide
/**
 * @author viktor
 */
object FromTramToTable {

  class PC(
    val istart: Int,
    val iapex: Int,
    val iend: Int,
    val rtPc: ArrayBuffer[RTPc])

  case class RTPc(rStart: Double, rApex: Double, rEnd: Double, area: Double)
  case class Datum(rt: Double, intensity: Double, id: String)
  case class Trace(id: String, rts: Seq[Double], intensities: Seq[Double])
  case class GroupedPeaks(id: String, peaks: ArrayBuffer[IdPeak])
  case class IdPeak(id: String, peaks: Seq[PC])
  case class Pair(pc: PC, id: String)
  case class PeakAndFrags(peak: PC, frags: ArrayBuffer[String])

  val traceIds = new ArrayBuffer[String]
  val ppm = 20 / 1e6
  var requests = new Queue[MasterRequest]
  var infinitePoller: ActorRef = _
  var printerActor: ActorRef = _
  var database = new ArrayBuffer[Trace]
  var databaseSize = 0
  var boundaryMin = 30
  var boundaryMax = 30
  var reps = 0
  var correlationLimit = 0.9
  var writeToFile = ""
  var results = new ArrayBuffer[Peak]
  var filter: Filter = _
  val comparer = new CompareToManual
  var detectAndGrouper: DetectAndGroup = _
  val filePrinter = new FilePrinter
  var group = new ArrayBuffer[GroupedPeaks]
  var thresholdOn = false
  val readRegression = new ReadRegression
  def setParams(min: Int, max: Int, detect: Int, filt: Int, rep: Int, ifTemplate: Int, corrLimit: Double, ifThresholdOn: Int) {
    boundaryMin = min
    boundaryMax = max
    correlationLimit = corrLimit
    detectAndGrouper = detect match {
      case 1 =>
        if (ifTemplate == 1) {

          SingleGroup(BruteForceComparer, TemplateDetect(correlationLimit))
        } else {

          SingleGroup(BruteForceComparer, SingleDerivativeDetect(boundaryMin, boundaryMax))
        }
      case 2 =>
        if (ifTemplate == 1) {

          SingleGroup(SortComparer, TemplateDetect(correlationLimit))
        } else {

          SingleGroup(SortComparer, SingleDerivativeDetect(boundaryMin, boundaryMax))
        }
      case 3 => MultiGroup(MultiDerivativeDetect(boundaryMin, boundaryMax))
    }
    filter = filt match {
      case 0 => NoFilter
      case 1 => SavitzkyGolay9
      case 2 => Daubchies2
      case 3 => Daubchies3
      case 4 => Daubchies8
      case 5 => Haar
      case 6 => ChainWavelet
      case 7 => MultiWavelet
    }

    reps = rep

    type MyType = Int => Boolean
    val equalToOne: MyType = (x: Int) => x == 1
    thresholdOn = equalToOne(ifThresholdOn)
  }
  def main(args: Array[String]) = {

    val tramlFile = new File(args(0))
    if (args.length >= 8) {
      setParams(args(1).toInt, args(2).toInt, args(3).toInt, args(4).toInt, args(5).toInt, args(6).toInt, args(7).toDouble, args(8).toInt)

    } else {
      filter = NoFilter
      detectAndGrouper = SingleGroup(SortComparer, TemplateDetect(correlationLimit))
    }

    println("reading traml")
    readTraml(tramlFile)
    val system = ActorSystem()
    //    val server = new InetSocketAddress("130.235.249.157", 33333)

    val server = new InetSocketAddress("localhost", 12345)
    infinitePoller = system.actorOf(Props[TramlPoller])
    printerActor = system.actorOf(Props[PrinterActor])
    println("connecting to server")
    val client = system.actorOf(MSDataProtocolActors.ClientInitiator.props(server, infinitePoller))
    system.awaitTermination
  }

  class PrinterActor extends Actor {
    var i = 0
    def receive = {
      case msg: Traces =>
        val data = new ArrayBuffer[Datum]
        for (pTrace <- msg.getPrecursorList) {
          val id = "prec %.4f".format(pTrace.getPrecursor.getLmz)
          database += Trace(id, pTrace.getTrace.getTimeList.asScala.map(_.doubleValue).toSeq, pTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
          println(id)
        }
        for (fTrace <- msg.getFragmentList) {

          val id = "frag %.4f -> %.4f".format(fTrace.getFragment.getPrecursor.getLmz, fTrace.getFragment.getFragment.getLmz)
          database += Trace(id, fTrace.getTrace.getTimeList.asScala.map(_.doubleValue).toSeq, fTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
          println(id)
        }
        if (databaseSize == database.size) {

          var featCalculator = new CalculateFeatures
          val time = System.currentTimeMillis()
          readRegression.readRegFile("coefFile")
          println("starting peak detection")
          var groupedFrags = new ArrayBuffer[Trace]
          var fragNames = new ArrayBuffer[String]
          var old = ""
          println(database.size)
          for (i <- 0 until database.size) {
            
            
            val ids = traceIds(i).split('/')
            val peptideId = ids(0) + "/" + ids(1).charAt(0)
            val fragId = ids(1).split('_')(1)
            if ((old != peptideId || database.size - 1 == i) && i != 0) {
              //calc features and do the rest
              var features = featCalculator.calcFeatures(groupedFrags, fragNames)(0)
              var min = Int.MaxValue
              var indexHolder = 0
              for (j <- 0 until readRegression.regressions.size) {
                val noOfPeaks = readRegression.estimatePeaks(features, readRegression.regressions(j))
                if (noOfPeaks < min && noOfPeaks > 0) {
                  min = noOfPeaks
                  indexHolder = j
                }
              }
              setParams(readRegression.regressions(indexHolder).min, readRegression.regressions(indexHolder).max, readRegression.regressions(indexHolder).detect, readRegression.regressions(indexHolder).filter, readRegression.regressions(indexHolder).rep, readRegression.regressions(indexHolder).ifTemp, readRegression.regressions(indexHolder).corrLim, readRegression.regressions(indexHolder).ifThresholdOn)

              results.appendAll(detectAndGrouper.doDetectionAndGrouping(groupedFrags, fragNames, filter, reps, thresholdOn))
              DetectionAndGrouping.clear()
              groupedFrags.clear()
              fragNames.clear()
            }
            groupedFrags += database(i)
            fragNames += traceIds(i)
            old = peptideId
          }
          //calc features
          //check regression
          //detect with least estimated peaks
          //results = detectAndGrouper.doDetectionAndGrouping(database, traceIds, filter, reps, thresholdOn)

          println("comparing with manualresult")

          var features = featCalculator.calcFeatures(database, traceIds)
          val manualResults = comparer.compareToCorrect(results, features, writeToFile)
          filePrinter.printToFile(boundaryMin, boundaryMax, detectAndGrouper, filter, correlationLimit, reps, manualResults, comparer.peptideHits, comparer.relevantResults, comparer.sensitivity, comparer.accuracy, writeToFile, comparer.fileName, thresholdOn, time)
          //          filePrinter.printFeatures(features, comparer.fileName)
          println("done")
          System.exit(0)
        }
    }
  }
  class TramlPoller extends Actor {

    import MSDataProtocolActors._

    def receive = {
      case msg: String =>
        println(msg)

      case MSDataProtocolConnected(remote, local) =>
        println("connected")
        if (requests.isEmpty) {
        } else {
          sender ! requests.dequeue
        }

      case MSDataReply(msg, nBytes, checkSum, timeTaken, remote) =>
        println("data recieved")
        printerActor ! msg.getTraces()
        if (requests.isEmpty) {
          println("closing")
        } else {
          sender ! requests.dequeue
        }

    }
  }

  def readTraml(f: File) = {
    val r = new XmlReader(new BufferedReader(new FileReader(f)))
    val traml = GhostTraML.fromFile(r)
    for (((pep, z), gts) <- traml.transitions.groupBy(gt => (gt.peptideRef, gt.q1z))) {
      for (gt <- gts) {
        traceIds += gt.id
        val req = GetTracesFor.newBuilder
        val prec = gt.q1
        val diff1 = prec * ppm
        val frag = gt.q3
        val diff2 = frag * ppm
        req.addFragment(
          FragmentBounds.newBuilder
            .setPrecursor(Bounds.newBuilder.setLmz(prec - diff1).setHmz(prec + diff1))
            .setFragment(Bounds.newBuilder.setLmz(frag - diff2).setHmz(frag + diff2)))
        requests += MasterRequest.newBuilder.setGetTracesFor(req).build
      }
    }
    databaseSize = requests.size
  }
}