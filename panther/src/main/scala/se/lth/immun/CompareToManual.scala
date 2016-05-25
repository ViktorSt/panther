package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import java.util.Scanner
import se.lth.immun.DataModel.FragmentArea
import java.io.File
import se.lth.immun.DataModel.Peak
import java.io.PrintWriter
import java.io.FileWriter

/**
 * @author viktor
 */
class CompareToManual {

  val resultReader = new Scanner(new File("result"))

  var relevantResults = 0
  var fragmentHits = 0
  var peptideHits = 0

  var sensitivity = 0.0
  var accuracy = 0.0
  var fileName = ""
  def compareToCorrect(results: ArrayBuffer[Peak], features: ArrayBuffer[PeptideFeatures], writeToFile: String): ArrayBuffer[Peak] = {
    fileName = resultReader.nextLine()
    var manualResults = new ArrayBuffer[Peak]
    var oldName = ""
    var flag = false
    var fragAndAreaList = new ArrayBuffer[FragmentArea]
    var oldApex = 0.0
    val areaPrinter = new java.io.PrintWriter(new File("areas"))
    while (resultReader.hasNext()) {
      val line = resultReader.nextLine()
      val data = line.split("\t")
      if (data(2) == " #N/A") {
      } else {
        fragAndAreaList += FragmentArea(data(1), data(3).toDouble)
        if (data(0) != oldName || !resultReader.hasNext()) {
          if (flag) {
            manualResults += new Peak(oldName,
              0,
              0,
              0,
              0,
              oldApex,
              0,
              fragAndAreaList,
              true)
          }

          flag = true
        }
        oldName = data(0)
        oldApex = data(2).toDouble
      }
    }
    flag = false
    val ifAssayFoundPrinter = new java.io.PrintWriter(new java.io.FileWriter("assayFound" + writeToFile))
    val testWriter = new java.io.PrintWriter(new java.io.FileWriter("validationTest"))
    ifAssayFoundPrinter.println(writeToFile)
    testWriter.println(manualResults.size)
    for (i <- 0 until manualResults.size) {

      var counter = 0
      for (j <- 0 until results.size) {
        if (manualResults(i).sequence == results(j).sequence) {
          counter += 1
          relevantResults += 1
          if (manualResults(i).isTheSame(results(j))) {
            if(flag==false){
            testWriter.println(manualResults(i).sequence + " " + manualResults(i).rtApex + " " + results(j).rtApex)
            }
            fragmentHits += 1
            flag = true
            results(j).truePeak = true
          }
        }
      }
      ifAssayFoundPrinter.println(manualResults(i).sequence + " " + flag + " " + counter)
      if (flag) {
        for (k <- 0 until features.size) {

          if (manualResults(i).sequence == features(k).peptideName) {
            features(k).peakFound = true
          }
        }

        peptideHits += 1
        flag = false
      }
    }
    ifAssayFoundPrinter.close()
    testWriter.close()
    if (peptideHits == 0 && fragmentHits == 0) {
    } else {
      sensitivity = peptideHits / manualResults.size.toDouble
      accuracy = fragmentHits / relevantResults.toDouble
    }
    return manualResults
  }

}