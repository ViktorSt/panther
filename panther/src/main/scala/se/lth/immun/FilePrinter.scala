package se.lth.immun

import java.util.Scanner
import java.io.File
import scala.collection.mutable.ArrayBuffer
import se.lth.immun.DataModel.Peak
import se.lth.immun.Filters.Filter
import se.lth.immun.DetectionAndGrouping._
import se.lth.immun.PeakFinder._
/**
 * @author viktor
 */
class FilePrinter {
  def printToFile(boundaryMin: Int, boundaryMax: Int, peakDetectionAndGrouping: DetectAndGroup, filter: Filter, correlationLimit: Double, reps: Int, manualResults: ArrayBuffer[Peak], peptideHits: Int, relevantResults: Int, sensitivity: Double, accuracy: Double, writeToFile: String, sample: String,thresholdOn: Boolean, startTime: Double) {
    val time = System.currentTimeMillis() - startTime
    val resultReader = new Scanner(new File("result"))
    val areaPrinter = new java.io.PrintWriter(new File("areas"))
    val plotWriter = new java.io.PrintWriter(new java.io.FileWriter(writeToFile + "plot", true))
    var filSiffra = 0
    if (sample == "napedro_L120227_008_SW") {
      filSiffra = 1
    } else if (sample == "napedro_L120302_005_SW") {
      filSiffra = 2
    } else if (sample == "napedro_L120420_004_SW") {
      filSiffra = 3
    }
    peakDetectionAndGrouping match {
      case SingleGroup(BruteForceComparer, TemplateDetect(_))            => plotWriter.println(filSiffra + "\t" + "n/a" + "\t" + "n/a" + "\t" + "Brute Force" + "\t" + "Template" + "\t" + correlationLimit + "\t" + filter + "\t" + sensitivity + "\t" + accuracy + "\t" + reps + "\t" + thresholdOn+"\t"+ time)

      case SingleGroup(BruteForceComparer, SingleDerivativeDetect(_, _)) => plotWriter.println(filSiffra + "\t" + boundaryMin + "\t" + boundaryMax + "\t" + "Brute Force" + "\t" + "Single Derivative" + "\t" + "n/a" + "\t" + filter + "\t" + sensitivity + "\t" + accuracy + "\t" + reps + "\t"+ thresholdOn+"\t" + time)
      case SingleGroup(SortComparer, TemplateDetect(_))                  => plotWriter.println(filSiffra + "\t" + "n/a" + "\t" + "n/a" + "\t" + "Sorted" + "\t" + "Template" + "\t" + correlationLimit + "\t" + filter + "\t" + sensitivity + "\t" + accuracy + "\t" + reps + "\t" + thresholdOn+"\t"+ time)
      case SingleGroup(SortComparer, SingleDerivativeDetect(_, _))       => plotWriter.println(filSiffra + "\t" + boundaryMin + "\t" + boundaryMax + "\t" + "Sorted" + "\t" + "Single Derivative" + "\t" + "n/a" + "\t" + filter + "\t" + sensitivity + "\t" + accuracy + "\t" + reps + "\t"+ thresholdOn+"\t" + time)
      case MultiGroup(MultiDerivativeDetect(_, _))                       => plotWriter.println(filSiffra + "\t" + boundaryMin + "\t" + boundaryMax + "\t" + "Multi" + "\t" + "Multi" + "\t" + "n/a" + "\t" + filter + "\t" + sensitivity + "\t" + accuracy + "\t" + reps + "\t"+ thresholdOn+"\t" + time)
    }
    plotWriter.close
  }
  def printFeatures(featureList: ArrayBuffer[PeptideFeatures], sample: String) {
    var filSiffra = 0
    if (sample == "napedro_L120227_008_SW") {
      filSiffra = 1
    } else if (sample == "napedro_L120302_005_SW") {
      filSiffra = 2
    } else if (sample == "napedro_L120420_004_SW") {
      filSiffra = 3
    }
    val featureWriter = new java.io.PrintWriter(new java.io.FileWriter("features" + filSiffra))
    for (i <- 0 until featureList.size) {
      featureWriter.println(featureList(i).peptideName + "\t" + featureList(i).average + "\t" + featureList(i).lowpassaverage + "\t" + featureList(i).fakelowpassaverage + "\t" + featureList(i).max + "\t" + featureList(i).median + "\t" + featureList(i).numberOfChanges + "\t" + featureList(i).variance)
    }
    featureWriter.close()
  }
  def printPeaksAndIfFound(manualResults: ArrayBuffer[Peak], results: ArrayBuffer[Peak]) {
    for (mr <- manualResults) {
      var flag = true
      for (i <- 0 until results.size) {
        if (mr.sequence == results(i).sequence) {
          println(results(i).sequence + " " + results(i).rtApex + results(i).truePeak)
          flag = false
        }
      }
      if (flag) {
        println(mr.sequence + " not found")
      }
    }
  }
}