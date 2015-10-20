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
import scala.collection.{ JavaConversions => JConversions };
import scala.collection.mutable.Queue
import se.lth.immun.traml.ghost.GhostTraML
import scala.collection.JavaConversions._
import collection.JavaConverters._

import akka.io.Tcp
import Tcp._
/**
 * @author viktor
 */
object FromTramToTable {

  class PC(
    val istart: Int,
    val iapex: Int,
    val iend: Int)
  def getDerivate(x: Seq[Double]): Seq[Double] =
    return x.zip(x.tail).map(tu => tu._2 - tu._1)
  case class ResultPc(rStart: Double, rApex: Double, rEnd: Double, area: Double)
  case class Datum(rt: Double, intensity: Double, id: String)
  case class Trace(id: String, rts: Seq[Double], intensities: Seq[Double])
  case class GroupedPeaks(id: String, peaks: ArrayBuffer[IdPeak])
  case class IdPeak(id: String, peaks: Seq[PC])
  case class PeakPair(peak: ResultPc, frags: ArrayBuffer[String])
  case class Result(assayId: String, peakPairs: ArrayBuffer[PeakPair], fragId: String)
  case class Pair(pc: PC, id: String)
  case class PeakAndFrags(peak: PC, frags: ArrayBuffer[String])
  case class ResultPAF(assayId: String, peakPairs: ArrayBuffer[PeakPair])

  val traceIds = new ArrayBuffer[String]
  val ppm = 20 / 1e6
  var sb = new Queue[MasterRequest]
  var infinitePoller: ActorRef = _
  var printerActor: ActorRef = _
  //  val output = new File("test")
  //  val p = new java.io.PrintWriter(output)
  var database = new ArrayBuffer[Trace]
  var peaks = new ArrayBuffer[IdPeak]
  var size = 0
  var group = new ArrayBuffer[GroupedPeaks]
  var results = new ArrayBuffer[Result]
  var resultsPAF = new ArrayBuffer[ResultPAF]
  var sample = ""
  var boundaryMin = 20
  var boundaryMax = 20
  var peakDetectionAndGrouping = 2
  var filter = 0
  var reps = 0

  val haar = Array(
    0.7071067811865475,
    0.7071067811865475)

  val gaar = Array(
    0.7071067811865475,
    -0.7071067811865475)

  val G8 = Array(
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
  val H8 = Array(
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
  val H3 = Array(
    math.sqrt(2) * (1 + math.sqrt(10) + math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (5 + math.sqrt(10) + 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (10 - 2 * math.sqrt(10) + 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (10 - 2 * math.sqrt(10) - 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (5 + math.sqrt(10) - 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (1 + math.sqrt(10) - math.sqrt(5 + 2 * math.sqrt(10))) / 32)
  val G3 = Array(H3(5), -H3(4), H3(3), -H3(2), H3(1), -H3(0))

  val H2 = Array(
    (1 + math.sqrt(3)) / (4 * math.sqrt(2)),
    (3 + math.sqrt(3)) / (4 * math.sqrt(2)),
    (3 - math.sqrt(3)) / (4 * math.sqrt(2)),
    (1 - math.sqrt(3)) / (4 * math.sqrt(2)))

  val G2 = Array(H2(3), -H2(2), H2(1), -H2(0))

  def main(args: Array[String]) = {
    if (args.length == 5) {
      boundaryMin = args(0).toInt
      boundaryMax = args(1).toInt
      peakDetectionAndGrouping = args(2).toInt
      filter = args(4).toInt
      reps = args(5).toInt
    }
    val temp = new File(args(0))
    println("reading traml")
    readTraml(temp)
    val scan = new Scanner(sb.toString());
    val system = ActorSystem()
    //    val server = new InetSocketAddress("130.235.249.157", 33333)

    val server = new InetSocketAddress("localhost", 12345)
    infinitePoller = system.actorOf(Props[TramlPoller])
    printerActor = system.actorOf(Props[PrinterActor])
    println("connecting to server")
    val client = system.actorOf(MSDataProtocolActors.ClientInitiator.props(server, infinitePoller))
    system.awaitTermination
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

  def findPCs(
    y: Seq[Double],
    dy: Seq[Double],
    ddy: Seq[Double]): Seq[PC] = {
    var ret = new ArrayBuffer[PC]

    var istart = 0
    var iapex = -1
    for (i <- 0 until dy.length - 1) {
      val starty = y(istart)
      if (math.signum(dy(i)) != math.signum(dy(i + 1))) {
        if (ddy(i) > 0) { // local min
          if (iapex >= 0 && boundaryMin >= y(i + 1) - starty) {
            //            println("valid min")
            var max = y(iapex)
            val endy = y(i + 1)
            if (boundaryMax <= max - math.max(starty, endy)) {
              //              println("found peak")
              ret += new PC(istart, iapex, i + 1)
            }

            istart = i + 1

            iapex = -1
          } else if (y(i + 1) < starty) {
            istart = i + 1
          }
        } else { // local max
          if (iapex == (-1)) {
            //            println("setting max, old apex -1")
            iapex = i + 1
          } else if (y(iapex) < y(i + 1)) {

            //            println("setting max, new apex higher")
            iapex = i + 1
          }

        }
      } else if (y(i) <= starty && i == istart + 1) { //concider if best approach
        istart = i
      }
    }

    return ret
  }

  class PrinterActor extends Actor {
    var i = 0
    def receive = {
      case msg: Traces =>
        val data = new ArrayBuffer[Datum]
        for (pTrace <- msg.getPrecursorList) {
          val id = "prec %.4f".format(pTrace.getPrecursor.getLmz)
          //                    val trace = pTrace.getTrace.getTimeList.zip(pTrace.getTrace.getIntensityList)
          database += Trace(id, pTrace.getTrace.getTimeList.asScala.map(_.doubleValue).toSeq, pTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
          //                    data ++= trace.map(t => Datum(t._1, t._2, id))
        }
        for (fTrace <- msg.getFragmentList) {

          val id = "frag %.4f -> %.4f".format(fTrace.getFragment.getPrecursor.getLmz, fTrace.getFragment.getFragment.getLmz)
          //                    val trace = fTrace.getTrace.getTimeList.zip(fTrace.getTrace.getIntensityList)
          //          if(i<4)
          //          println(fTrace.getTrace.getTimeList.size + " " + fTrace.getTrace.getIntensityList.size + " " + id + " " + traceIds(i))
          //          i+=1
          database += Trace(id, fTrace.getTrace.getTimeList.asScala.map(_.doubleValue).toSeq, fTrace.getTrace.getIntensityList.asScala.map(_.doubleValue).toSeq)
          //                    data ++= trace.map(t => Datum(t._1, t._2, id))
        }
        //                for (d <- data) {
        //                  val temp = d.id + '\t' + d.rt + '\t' + d.intensity
        //                  p.println(temp)
        //                }
        if (size == database.size) {

          val time = System.currentTimeMillis()
          println("starting peak detection")
          if (peakDetectionAndGrouping == 1) {
            groupData()
            println("comparing peaks")
            comparePeaks()
          } else if (peakDetectionAndGrouping == 2) {
            groupData()
            println("comparing peaks")
            newCompare()
          } else if (peakDetectionAndGrouping == 3) {
            newGrouping()
          }
          //          printResult()
          println("comparing with manualresult")
          compareToCorrect((System.currentTimeMillis() - time))
        }
    }
  }
  def compareToCorrect(time: Double) {
    var i = 0
    var j = 0
    var count = 0

    val scan = new Scanner(new File("result"))
    val p = new java.io.PrintWriter(new File("comparison"))
    sample = scan.nextLine()
    while (scan.hasNext()) {
      val temp = scan.nextLine()
      val temp2 = temp.split("\t")
      var flag = false
      if (temp2(2) == " #N/A") {
      } else {
        j += 1
        for (t <- results) {
          if (t.assayId == temp2(0) && t.fragId == temp2(1)) {
            for (pp <- t.peakPairs) {
              count += 1
              if (pp.peak.rApex < temp2(2).toDouble * (1.001) && pp.peak.rApex * (1.001) > temp2(2).toDouble) {
                flag = true
                p.print("Found peak for " + temp2(0) + " " + temp2(1) + " " + pp.peak.rStart + " " + pp.peak.rApex + " " + pp.peak.rEnd + " " + temp2(3) + " " + pp.peak.area + " ")
                for (f <- pp.frags) {
                  p.print(f + " ")
                }
                p.println()
              }
            }
          }
        }
        for (t <- resultsPAF) {
          if (t.assayId == temp2(0)) {
            for (z <- 0 until t.peakPairs.length) {

              count += 1
              if (t.peakPairs(z).frags.contains(temp2(1))) {
                if (t.peakPairs(z).peak.rApex < temp2(2).toDouble * (1.01) && t.peakPairs(z).peak.rApex * (1.01) > temp2(2).toDouble) {
                  flag = true
                  p.print("Found peak for " + temp2(0) + " " + temp2(1) + " " + t.peakPairs(z).peak.rStart + " " + t.peakPairs(z).peak.rApex + " " + t.peakPairs(z).peak.rEnd + " " + temp2(3) + " " + t.peakPairs(z).peak.area + " ")
                  for (f <- t.peakPairs(z).frags) {
                    p.print(f + " ")
                  }
                  p.println()
                }
              }
            }
          }
        }
        p.println()
        if (!flag) {
          i += 1
          p.println("Peak not found for " + temp2(0) + " " + temp2(1))
        }
      }
    }

    p.close()

    val z = new java.io.PrintWriter(new java.io.FileWriter("data", true))
    //Skriv till en fil istället
    z.println(sample)
    z.println("boundaryMin : " + boundaryMin + " boundaryMax : " + boundaryMax)
    z.println("peakDetectionAndGrouping method: " + peakDetectionAndGrouping)
    z.println("Number of peaks not found: " + i + " of " + j + " possible" + " and " + count + " candidates")
    val truePos = (j - i) / j.toDouble
    val falsePos = (count - (j - i)) / count.toDouble
    z.println("Finding rate: " + truePos)
    z.println("False positive rate: " + falsePos)
    z.println("Time elapsed: " + time + "ms")
    z.println()
    z.close
    println("DONE")
  }
  def printResult() {
    for (i <- 0 until results.size) {
      println(results(i).assayId)
      for (j <- 0 until results(i).peakPairs.size) {
        println(j + "\n" + results(i).peakPairs(j).peak.rStart + " " + results(i).peakPairs(j).peak.rApex + " " + results(i).peakPairs(j).peak.rEnd)
        for (k <- 0 until results(i).peakPairs(j).frags.size) {
          print(results(i).peakPairs(j).frags(k) + '\t')
        }
        println()
      }
    }
  }

  def comparePeaks() {
    for (i <- 0 until group.size) {
      val g1 = group(i)
      for (j <- 0 until g1.peaks.size) {
        val p = g1.peaks(j)
        for (k <- 0 until p.peaks.size) {
          val pc = p.peaks(k)
          var temp = new ArrayBuffer[String]
          for (l <- 0 until g1.peaks.size) {

            var flag = false
            if (j != l) {
              flag = checkMatch(pc, g1.peaks(l).peaks)
            }
            if (flag) {
              temp += g1.peaks(l).id

            }
          }
          if (!temp.isEmpty) {
            var newPc = ResultPc(0, 0, 0, 0)
            for (l <- 0 until database.size) {
              val temp1 = traceIds(l).split('/')
              val temp2 = temp1(0) + "/" + temp1(1).charAt(0)
              if (temp2 == g1.id && p.id == temp1(1).split('_')(1)) {
                var sum = 0.0
                for (z <- pc.istart to pc.iend) {
                  if (z == 0) {

                    sum += database(l).intensities(z) * (database(l).rts(z)) * 60
                  } else {
                    sum += database(l).intensities(z) * (database(l).rts(z) - database(l).rts(z - 1)) * 60
                  }
                }

                //              println(temp2 + temp1(1).split('_')(1) + " " + database(l).rts.size + " " + pc.iend)
                //                                      println(pc.istart + " " + pc.iapex + " " + pc.iend)
                //                                      println(database(l).rts.size)
                //                      println("||||||||||||||")
                //                      println(database(l).rts(pc.istart) + " " + database(l).rts(pc.iapex) + " " + database(l).rts(pc.iend) +"\n\n\n")
                newPc = ResultPc(database(l).rts(pc.istart), database(l).rts(pc.iapex), database(l).rts(pc.iend), sum)
              }
            }
            //            if (temp.size == g1.peaks.size - 1) {
            var contains = false
            for (l <- 0 until results.size) {
              if (results(l).assayId == g1.id && results(l).fragId == p.id) {
                contains = true
                results(l).peakPairs += PeakPair(newPc, p.id +: temp)
              }
            }
            if (!contains) {
              results += Result(g1.id, new ArrayBuffer += PeakPair(newPc, p.id +: temp), p.id)
            }
            //            }
          }
        }
      }
    }
  }
  def newCompare() {
    for (i <- 0 until group.size) {
      val g1 = group(i)
      val temp = groupSorted(sortGroup(g1.peaks))
      if (!temp.isEmpty) {
        for (j <- 0 until temp.size) {
          val pc = temp(j).peak
          var newPc = ResultPc(0, 0, 0, 0)
          for (l <- 0 until database.size) {
            val temp1 = traceIds(l).split('/')
            val temp2 = temp1(0) + "/" + temp1(1).charAt(0)
            if (temp2 == g1.id && temp(j).frags(0) == temp1(1).split('_')(1)) {
              var sum = 0.0
              for (z <- pc.istart to pc.iend) {
                if (z == 0) {

                  sum += database(l).intensities(z) * (database(l).rts(z)) * 60
                } else {
                  sum += database(l).intensities(z) * (database(l).rts(z) - database(l).rts(z - 1)) * 60
                }
              }

              //              println(temp2 + temp1(1).split('_')(1) + " " + database(l).rts.size + " " + pc.iend)
              //                                      println(pc.istart + " " + pc.iapex + " " + pc.iend)
              //                                      println(database(l).rts.size)
              //                      println("||||||||||||||")
              //                      println(database(l).rts(pc.istart) + " " + database(l).rts(pc.iapex) + " " + database(l).rts(pc.iend) +"\n\n\n")
              newPc = ResultPc(database(l).rts(pc.istart), database(l).rts(pc.iapex), database(l).rts(pc.iend), sum)
            }
          }
          //            if (temp.size == g1.peaks.size - 1) {
          var contains = false
          for (l <- 0 until resultsPAF.size) {
            if (resultsPAF(l).assayId == g1.id) {
              contains = true
              resultsPAF(l).peakPairs += PeakPair(newPc, temp(j).frags)
            }
          }
          if (!contains) {
            resultsPAF += new ResultPAF(g1.id, new ArrayBuffer += PeakPair(newPc, temp(j).frags))
          }

          //            }
        }
      }
    }
  }
  //Används ej?
  def sortGroup(peaks: ArrayBuffer[IdPeak]): ArrayBuffer[Pair] = {
    var i = 0
    val temp = new ArrayBuffer[Pair]
    while (i < peaks.size) {
      var j = 0
      while (j < peaks(i).peaks.size) {
        if (temp.isEmpty) {
          temp += Pair(peaks(i).peaks(j), peaks(i).id)
        }
        for (k <- 0 until temp.size) {
          if (temp(k).pc.istart > peaks(i).peaks(j).istart) {
            temp.insert(k, Pair(peaks(i).peaks(j), peaks(i).id))
          } else if (k == temp.size - 1 && temp(k).pc.istart < peaks(i).peaks(j).istart) {
            temp += (Pair(peaks(i).peaks(j), peaks(i).id))
          }
        }
        j += 1
      }
      i += 1
    }

    return temp
  }

  //Används ej heller
  def groupSorted(peaks: ArrayBuffer[Pair]): ArrayBuffer[PeakAndFrags] = {
    var i = 0
    var temp = new ArrayBuffer[PeakAndFrags]
    while (i < peaks.size) {
      var flag = true
      var j = i
      var max = peaks(i).pc.iend
      var frags = new ArrayBuffer[String]
      while (flag) {

        j += 1
        if (j >= peaks.size) {
          flag = false
        } else if (peaks(i).pc.iend >= peaks(j).pc.istart) {
          if (!frags.contains(peaks(j).id)) {
            frags += peaks(j).id
          }
          //add to temp
          if (peaks(j).pc.iend > max) {
            max = peaks(j).pc.iend
          }
        } else {
          flag = false

        }
      }
      if (frags.size > 1) {
        temp += new PeakAndFrags(new PC(peaks(i).pc.istart, peaks(i).pc.iapex, max), frags)

      }
      i = j
      //if temp tillräckligt stor add to result
    }
    return temp
  }

  def chainWavelet(x: Seq[Double], levels: Int): Seq[Double] = {
    var yh1 = new Array[Double](x.size / 2)
    var yl1 = new Array[Double](x.size / 2)
    val threshhold = 0.1
    var temp1 = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform
    {
      yh1(i) = x((i * 2) % x.size) * gaar(1) + x((i * 2 + 1) % x.size) * gaar(0)
      if (yh1(i) < 10) {
        yh1(i) = 0
      }
      yl1(i) = x((i * 2) % x.size) * haar(1) + x((i * 2 + 1) % x.size) * haar(0)
    }

    var yh2 = new Array[Double](yl1.size / 2)
    var yl2 = new Array[Double](yl1.size / 2)

    var temp2 = new Array[Double](yl1.size)
    for (i <- 1 until yl1.size / 2) //forward Transform
    {
      yh2(i) = yl1((i * 2) % yl1.size) * G2(3) + yl1((i * 2 + 1) % yl1.size) * G2(2) + yl1((i * 2 + 2) % yl1.size) * G2(1) + yl1((i * 2 + 3) % yl1.size) * G2(0)
      if (yh2(i) < 10) {
        yh2(i) = 0
      }
      yl2(i) = yl1((i * 2) % yl1.size) * H2(3) + yl1((i * 2 + 1) % yl1.size) * H2(2) + yl1((i * 2 + 2) % yl1.size) * H2(1) + yl1((i * 2 + 3) % yl1.size) * H2(0)
    }

    if (levels > 0) {
      chainWavelet(yl2, levels - 1)
    }
    for (i <- 0 until yl1.size / 2) //inverse Transform
    {
      temp2((i * 2) % x.size) = yl2(((i - 1) % yl2.size + yl2.size) % yl2.size) * H2(1) + yl2((i) % yl2.size) * H2(3) + yh2((((yh2.size - 1 + i) % yh2.size) + yh2.size) % yh2.size) * G2(1) + yh2((yh2.size + i) % yh2.size) * G2(3)
      temp2((i * 2 + 1) % x.size) = yl2(((yl2.size - 1 + i) % yl2.size + yl2.size) % yl2.size) * H2(0) + yl2((yl2.size + i) % yl2.size) * H2(2) + yh2(((yh2.size - 1 + i) % yh2.size + yh2.size) % yh2.size) * G2(0) + yh2((yh2.size + i) % yh2.size) * G2(2)

    }

    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp1((i * 2) % x.size) = temp2(((i) % temp2.size + temp2.size) % temp2.size) * haar(1) + yh1((((yh1.size - 1 + i) % yh1.size) + yh1.size) % yh1.size) * gaar(1)
      temp1((i * 2 + 1) % x.size) = temp2(((temp2.size + i) % temp2.size + temp2.size) % temp2.size) * haar(0) + yh1(((yh1.size + i) % yh1.size + yh1.size) % yh1.size) * gaar(0)

    }

    return temp1.toSeq map (d => math.max(0, d))
  }

  def multiWavelet(x: Seq[Double], levels: Int): Seq[Double] = {
    var temp = daubchies2(x, 0)
    temp = daubchies8(temp, 0)
    if (levels > 0) {
      temp = multiWavelet(temp, levels - 1)
    }

    return temp
  }

  def haar(x: Seq[Double], levels: Int): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    val threshhold = 0.1
    var temp = new Array[Double](x.size)
    yh(0) = 0
    yl(0) = x(0)
    for (i <- 1 until x.size / 2) //forward Transform
    {
      yh(i) = x((i * 2) % x.size) * gaar(1) + x((i * 2 + 1) % x.size) * gaar(0)
      if (yh(i) < 10) {
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * haar(1) + x((i * 2 + 1) % x.size) * haar(0)
    }
    if (levels > 0) {
      yl = haar(yl.toSeq, levels - 1).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i) % yl.size + yl.size) % yl.size) * haar(1) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * gaar(1)
      temp((i * 2 + 1) % x.size) = yl(((yl.size + i) % yl.size + yl.size) % yl.size) * haar(0) + yh(((yh.size + i) % yh.size + yh.size) % yh.size) * gaar(0)

    }

    return temp.toSeq map (d => math.max(0, d))
  }

  def daubchies2(x: Seq[Double], levels: Int): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    val threshhold = 0.1
    var temp = new Array[Double](x.size)
    yh(0) = 0
    yl(0) = x(0)
    for (i <- 1 until x.size / 2) //forward Transform
    {
      yh(i) = x((i * 2) % x.size) * G2(3) + x((i * 2 + 1) % x.size) * G2(2) + x((i * 2 + 2) % x.size) * G2(1) + x((i * 2 + 3) % x.size) * G2(0)
      if (yh(i) < 10) {
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * H2(3) + x((i * 2 + 1) % x.size) * H2(2) + x((i * 2 + 2) % x.size) * H2(1) + x((i * 2 + 3) % x.size) * H2(0)
    }
    if (levels > 0) {
      yl = daubchies2(yl.toSeq, levels - 1).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 1) % yl.size + yl.size) % yl.size) * H2(1) + yl((i) % yl.size) * H2(3) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G2(1) + yh((yh.size + i) % yh.size) * G2(3)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H2(0) + yl((yl.size + i) % yl.size) * H2(2) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G2(0) + yh((yh.size + i) % yh.size) * G2(2)

    }

    return temp.toSeq map (d => math.max(0, d))
  }

  def daubchies3(x: Seq[Double], levels: Int): Seq[Double] = {
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
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform
    {

      yh(i) = x((i * 2) % x.size) * G8(15) + x((i * 2 + 1) % x.size) * G8(14) + x((i * 2 + 2) % x.size) * G8(13) + x((i * 2 + 3) % x.size) * G8(12) + x((i * 2 + 4) % x.size) * G8(11) + x((i * 2 + 5) % x.size) * G8(10) + x((i * 2 + 6) % x.size) * G8(9) + x((i * 2 + 7) % x.size) * G8(8) + x((i * 2 + 8) % x.size) * G8(7) + x((i * 2 + 9) % x.size) * G8(6) + x((i * 2 + 10) % x.size) * G8(5) + x((i * 2 + 11) % x.size) * G8(4) + x((i * 2 + 12) % x.size) * G8(3) + x((i * 2 + 13) % x.size) * G8(2) + x((i * 2 + 14) % x.size) * G8(1) + x((i * 2 + 15) % x.size) * G8(0)
      if (yh(i) < 10) {
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

  def newGrouping() {
    var old = ""
    var j = -1
    var tempGroup = new ArrayBuffer[Seq[Double]]
    for (i <- 0 until database.size) {
      val temp1 = traceIds(i).split('/')
      val temp2 = temp1(0) + "/" + temp1(1).charAt(0)

      var signalOrg = database(i).intensities
      if (filter == 1) {
        signalOrg = savitzkyGolay9(signalOrg)
      } else if (filter == 2) {

        signalOrg = daubchies2(signalOrg, reps)
      } else if (filter == 3) {

        signalOrg = daubchies3(signalOrg, reps)
      } else if (filter == 4) {

        signalOrg = daubchies8(signalOrg, reps)
      } else if (filter == 5) {

        signalOrg = haar(signalOrg, reps)
      } else if (filter == 6) {

        signalOrg = chainWavelet(signalOrg, reps)
      } else if (filter == 7) {

        signalOrg = multiWavelet(signalOrg, reps)
      }
      if (temp2 == old && i != database.size - 1) {
        tempGroup += signalOrg
      } else {
        if (!tempGroup.isEmpty) {
          var testPeaks = findPeaks(tempGroup, i)
          tempGroup = new ArrayBuffer[Seq[Double]]
          resultsPAF += new ResultPAF(old, new ArrayBuffer[PeakPair])
          for (t <- testPeaks) {
            resultsPAF(resultsPAF.length - 1).peakPairs += new PeakPair(new ResultPc(database(i - 1).rts(t.peak.istart), database(i - 1).rts(t.peak.iapex), database(i - 1).rts(t.peak.iend), 0), t.frags)

          }
        }

        tempGroup += signalOrg

        j += 1
        old = temp2
      }

    }
  }
  //Return PeakAndFrags istället
  def findPeaks(peakLists: ArrayBuffer[Seq[Double]], index: Int): Seq[PeakAndFrags] = {
    var listOfFrags = new ArrayBuffer[String]
    for (i <- peakLists.size to 1 by -1) {
      val temp1 = traceIds(index - i).split('/')
      listOfFrags += temp1(1).split('_')(1)
    }

    var ys = new ArrayBuffer[Seq[Double]]
    var dys = new ArrayBuffer[Seq[Double]]
    var ddys = new ArrayBuffer[Seq[Double]]
    var ret = new ArrayBuffer[PeakAndFrags]
    var minflag = false
    var maxflag = false

    var flags = Array(false, false, false, false)
    for (pl <- peakLists) {
      ys += pl
      dys += getDerivate(pl)
      ddys += getDerivate(getDerivate(pl))
    }
    var istart = 0
    var iapex = -1
    var starty = ys(0)(istart)
    var endy = -1.0
    var max = -1.0
    var old = true
    for (i <- 0 until dys(0).length - 1) {
      var j = -1
      var k = 0
      var z = -1
      var l = 0
      var tempList = new ArrayBuffer[String]
      for (t <- 0 until dys.length) {
        j += 1
        if (i == 0) {
          if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
            if (ddys(t)(i + 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        } else if (i == dys(0).length - 2) {
          if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
            if (ddys(t)(i - 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        } else {
          if (math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
            if (ddys(t)(i + 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
            if (ddys(t)(i - 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        }

      }

      if (l > 1) { // local max
        var temp = ys(0)(i + 1)
        for (y <- ys) {
          if (y(i + 1) > temp) {
            temp = y(i + 1)
          }
        }
        if (iapex == -1) {
          old = false
          iapex = i + 1
          max = temp
        } else if (temp > max) {
          old = false
          iapex = i + 1
          max = temp
        }
        for (t <- 0 until dys.length) {
          if (i == 0) {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
              flags(t) = true
            }
          } else if (i == dys(0).length - 2) {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
              flags(t) = true
            }
          } else {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2)) || math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
              flags(t) = true
            }
          }
        }

      } else if (k > 1) { // local min

        for (y <- ys) {
          if (boundaryMin >= y(i + 1) - starty) {
            minflag = true
          }
        }
        if (iapex >= 0 && minflag) {
          var temp = ys(0)(i + 1)
          for (y <- ys) {
            if (y(i + 1) < temp) {
              temp = y(i + 1)
            }
          }
          endy = temp
          if (boundaryMax <= max - math.max(starty, endy)) {
            var temp = new ArrayBuffer[String]
            for (t <- 0 until ddys.length) {
              if (flags(t)) {
                temp += listOfFrags(t)
              }
            }
            ret += new PeakAndFrags(new PC(istart, iapex, i + 1), temp)
          }

          flags = Array(false, false, false, false)
          iapex = -1

          istart = i + 1
          old = false
        }

        minflag = false

      } else if (ys(0)(i) < ys(0)(istart) && i == istart + 1) {
        istart = i
        old = true
      } else {
        old = true
      }
    }
    return ret
  }
  def groupData() {
    var old = ""
    var j = -1
    for (i <- 0 until database.size) {
      var contains = false
      val temp1 = traceIds(i).split('/')
      val temp2 = temp1(0) + "/" + temp1(1).charAt(0)
      var signalOrg = database(i).intensities
      if (filter == 1) {
        signalOrg = savitzkyGolay9(signalOrg)
      } else if (filter == 2) {

        signalOrg = daubchies2(signalOrg, reps)
      } else if (filter == 3) {

        signalOrg = daubchies3(signalOrg, reps)
      } else if (filter == 4) {

        signalOrg = daubchies8(signalOrg, reps)
      } else if (filter == 5) {

        signalOrg = haar(signalOrg, reps)
      } else if (filter == 6) {

        signalOrg = chainWavelet(signalOrg, reps)
      } else if (filter == 7) {

        signalOrg = multiWavelet(signalOrg, reps)
      }
      val temp = getDerivate(signalOrg)
      val seq = findPCs(signalOrg, temp, getDerivate(temp))
      if (temp2 == old) {
        group(j).peaks += IdPeak(temp1(1).split('_')(1), seq)
      } else {
        j += 1
        group += GroupedPeaks(temp2, new ArrayBuffer += IdPeak(temp1(1).split('_')(1), seq))
        old = temp2
      }

    }
  }
  def checkMatch(peak: PC, peaks: Seq[PC]): Boolean = {
    for (p <- peaks) {
      //      if (peak.istart < p.iapex && p.iapex < peak.iend) {
      if ((peak.istart <= p.istart && peak.iend >= p.istart) || (peak.istart <= p.iend && peak.iend >= p.iend) || (peak.istart >= p.istart && peak.iend <= p.iend) || (peak.istart <= p.istart && peak.iend >= p.iend)) {
        return true
      } else if (peak.iend < p.istart) {
        return false
      }
    }
    return false
  }
  class TramlPoller extends Actor {

    import MSDataProtocolActors._

    def receive = {
      case msg: String =>
        println(msg)

      case MSDataProtocolConnected(remote, local) =>
        println("connected")
        if (sb.isEmpty) {
        } else {
          sender ! sb.dequeue
        }

      case MSDataReply(msg, nBytes, checkSum, timeTaken, remote) =>
        println("data recieved")
        printerActor ! msg.getTraces()
        if (sb.isEmpty) {
          println("closing")
        } else {
          sender ! sb.dequeue
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
        sb += MasterRequest.newBuilder.setGetTracesFor(req).build
      }
    }
    size = sb.size
  }
}