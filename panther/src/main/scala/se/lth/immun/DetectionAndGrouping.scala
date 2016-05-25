package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import se.lth.immun.FromTramToTable.Trace
import se.lth.immun.FromTramToTable.PC
import se.lth.immun.FromTramToTable.IdPeak
import se.lth.immun.FromTramToTable.GroupedPeaks
import se.lth.immun.DataModel.Peak
import se.lth.immun.DataModel.FragmentArea
import se.lth.immun.FromTramToTable.Pair
import se.lth.immun.FromTramToTable.PeakAndFrags
import se.lth.immun.FromTramToTable.RTPc
import se.lth.immun.Filters.Filter
import se.lth.immun.PeakFinder._

/**
 * @author viktor
 */
object DetectionAndGrouping {
  trait DetectAndGroup{
  def doDetectionAndGrouping(database: ArrayBuffer[Trace], traceIds: ArrayBuffer[String], filter: Filter, reps: Int, thresholdOn: Boolean):ArrayBuffer[Peak]    
  
  }
  case class SingleGroup(compareMethod: CompareMethod, peakFinder: SinglePeakFinder) extends DetectAndGroup{
    def doDetectionAndGrouping(database: ArrayBuffer[Trace], traceIds: ArrayBuffer[String], filter: Filter, reps: Int, thresholdOn: Boolean):ArrayBuffer[Peak] = 
    {
       var old = ""
    var j = -1
    for (i <- 0 until database.size) {
      val ids = traceIds(i).split('/')
      val peptideId = ids(0) + "/" + ids(1).charAt(0)
      val fragId=ids(1).split('_')(1)
      var signal = filter.applyFilter(database(i).intensities, reps, thresholdOn)
      var seq = Seq[PC]()
        seq =peakFinder.findPc(signal, database(i).rts)
      if (peptideId == old) {
        group(j).peaks += IdPeak(fragId, seq)
      } else {
        j += 1
        group += GroupedPeaks(peptideId, new ArrayBuffer += IdPeak(fragId, seq))
        old = peptideId
      }
    }
    compareMethod.compare()
    
    return results
  }
  }
  case class MultiGroup(peakFinder: MultiDerivativeDetect) extends DetectAndGroup{
    def doDetectionAndGrouping(database: ArrayBuffer[Trace], traceIds: ArrayBuffer[String], filter: Filter, reps: Int, thresholdOn: Boolean):ArrayBuffer[Peak]=
      {
    var oldId = ""
    var j = -1
    var testmatrix = Array.fill(9)(0.0)
    var ys = new ArrayBuffer[Seq[Double]]
    for (i <- 0 until database.size) {
      val ids = traceIds(i).split('/')

      val peptideId = ids(0) + "/" + ids(1).charAt(0)
      
      var signal = filter.applyFilter(database(i).intensities, reps, thresholdOn)
      if (peptideId == oldId && i != database.size - 1) {
        ys += signal
      } else {
        if (!ys.isEmpty) {

          var listOfFrags = new ArrayBuffer[String]
          for (z <- ys.size to 1 by -1) {
            val id = traceIds(i - z).split('/')
            listOfFrags += id(1).split('_')(1)
          }
          var testPeaks = peakFinder.findPc(ys, listOfFrags, database(i).rts)
          ys = new ArrayBuffer[Seq[Double]]
          for (t <- testPeaks) {
            var fragmentAreas = new ArrayBuffer[FragmentArea]
            for (i <- 0 until t.frags.size) {
              fragmentAreas += FragmentArea(t.frags(i), t.peak.rtPc(i).area)
            }
            results += new Peak(oldId,
              t.peak.istart,
              t.peak.iapex,
              t.peak.iend,
              t.peak.rtPc(0).rStart,
              t.peak.rtPc(0).rApex,
              t.peak.rtPc(0).rEnd,
              fragmentAreas,
              false)
          }
        }

        ys += signal

        j += 1
        oldId = peptideId
      }

    }
    return results
  }

  }
  trait CompareMethod{
    def compare()
  }
  case object BruteForceComparer extends CompareMethod{
    def compare()=bruteForceCompare()
  }
  
  case object SortComparer extends CompareMethod{
    def compare()=sortCompare()
  }
  def clear(){
    results.clear()
    group.clear()
  }
  var singlePeakfinder :SinglePeakFinder=_
  var results=new ArrayBuffer[Peak]
  var group = new ArrayBuffer[GroupedPeaks]
  def bruteForceCompare() {
    for (i <- 0 until group.size) {
      val g1 = group(i)
      for (j <- 0 until g1.peaks.size) {
        val p = g1.peaks(j)
        for (k <- 0 until p.peaks.size) {
          val pc = p.peaks(k)
          var frags = new ArrayBuffer[String]
          for (l <- 0 until g1.peaks.size) {

            var flag = false
            if (j != l) {
              flag = checkMatch(pc, g1.peaks(l).peaks)
            }
            if (flag) {
              frags += g1.peaks(l).id

            }
          }
          if (!frags.isEmpty) {

            var fragArea = new ArrayBuffer[FragmentArea]
            fragArea += FragmentArea(p.id, pc.rtPc(0).area)
            for (z <- 1 until frags.size) {
              fragArea += FragmentArea(frags(z), 0)
            }
            results += new Peak(g1.id,
              pc.istart,
              pc.iapex,
              pc.iend,
              pc.rtPc(0).rStart,
              pc.rtPc(0).rApex,
              pc.rtPc(0).rEnd,
              fragArea,
              false)
          }
        }
      }
    }
  }
  
  def checkMatch(peak: PC, peaks: Seq[PC]): Boolean = {
    for (p <- peaks) {
      if ((peak.istart <= p.istart && peak.iend >= p.istart) || (peak.istart <= p.iend && peak.iend >= p.iend) || (peak.istart >= p.istart && peak.iend <= p.iend) || (peak.istart <= p.istart && peak.iend >= p.iend)) {
        return true
      } else if (peak.iend < p.istart) {
        return false
      }
    }
    return false
  }
  
  def sortCompare() {

    def mergeSortedPCs(pcs1: IdPeak, pcs2: ArrayBuffer[Pair]): ArrayBuffer[Pair] = {
      var i1 = 0
      var i2 = 0

      val res = new ArrayBuffer[Pair]
      while (i1 < pcs1.peaks.length && i2 < pcs2.length) {
        val p1 = pcs1.peaks(i1)
        val p2 = pcs2(i2)
        if (p1.istart < p2.pc.istart) {
          res += new Pair(p1, pcs1.id)
          i1 += 1
        } else {
          res += p2
          i2 += 1
        }
      }

      while (i1 < pcs1.peaks.length) {
        res += new Pair(pcs1.peaks(i1), pcs1.id)
        i1 += 1
      }

      while (i2 < pcs2.length) {
        res += pcs2(i2)
        i2 += 1
      }

      res
    }

    for (i <- 0 until group.size) {
      val g1 = group(i)

      val sortedG1 = g1.peaks.map(idPeak => new IdPeak(idPeak.id, idPeak.peaks.sortBy(_.istart)))

      var sortedPCs = new ArrayBuffer[Pair]
      for (j <- 0 until sortedG1.size) {
        sortedPCs = mergeSortedPCs(sortedG1(j), sortedPCs)
      }
      val tempResult = groupSorted(sortedPCs)
      if (!tempResult.isEmpty) {
        for (j <- 0 until tempResult.size) {
          val pc = tempResult(j).peak
          var fragArea = new ArrayBuffer[FragmentArea]
          fragArea += FragmentArea(tempResult(j).frags(0), pc.rtPc(0).area)
          for (z <- 1 until tempResult(j).frags.size) {
            fragArea += FragmentArea(tempResult(j).frags(z), 0)
          }
          results += new Peak(g1.id,
            pc.istart,
            pc.iapex,
            pc.iend,
            pc.rtPc(0).rStart,
            pc.rtPc(0).rApex,
            pc.rtPc(0).rEnd,
            fragArea,
            false)

        }
      }
    }
  }
  
  
  
  def groupSorted(peaks: ArrayBuffer[Pair]): ArrayBuffer[PeakAndFrags] = {
    var i = 0
    var ret = new ArrayBuffer[PeakAndFrags]
    while (i < peaks.size) {
      var flag = true
      var j = i
      var iMax = peaks(i).pc.iend
      var rMax = peaks(i).pc.rtPc(0).rEnd
      var frags = new ArrayBuffer[String] += peaks(i).id
      while (flag) {

        j += 1
        if (j >= peaks.size) {
          flag = false
        } else if (peaks(i).pc.iend >= peaks(j).pc.istart) {
          if (!frags.contains(peaks(j).id)) {
            frags += peaks(j).id
          }
          if (peaks(j).pc.iend > iMax) {
            iMax = peaks(j).pc.iend
            rMax = peaks(j).pc.rtPc(0).rEnd
          }
        } else {
          flag = false

        }
      }
      if (frags.size > 1) {
        ret += new PeakAndFrags(
          new PC(
            peaks(i).pc.istart,
            peaks(i).pc.iapex,
            iMax,
            ArrayBuffer(
              new RTPc(
                peaks(i).pc.rtPc(0).rStart,
                peaks(i).pc.rtPc(0).rApex,
                rMax,
                peaks(i).pc.rtPc(0).area))),
          frags)

      }
      i = j
    }
    return ret
  }
}