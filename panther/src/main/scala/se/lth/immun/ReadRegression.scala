package se.lth.immun

import java.util.Scanner
import java.io.File
import scala.collection.mutable.ArrayBuffer

/**
 * @author viktor
 */
class ReadRegression {
  val regressions = new ArrayBuffer[RegressionParams]

  def readRegFile(name: String) {
    val fileReader = new Scanner(new File(name))
    while (fileReader.hasNext()) {
      val line = fileReader.nextLine()

      val splitData = line.split("\\)")
      val min = splitData(0).charAt(splitData(0).length - 1)-'0'
      val max = splitData(1).charAt(splitData(1).length - 1)-'0'
      val detect = splitData(2).charAt(splitData(2).length - 1)-'0'
      val filter = splitData(3).charAt(splitData(3).length - 1)-'0'
      val rep = splitData(4).charAt(splitData(4).length - 1)-'0'
      val ifTemplate = splitData(5).charAt(splitData(5).length - 1)-'0'
      val corrLim = (splitData(6).split("\\(")(1)).toDouble
      val regParams = splitData(7).drop(2).split(" ").map { x => x.toDouble }

      regressions += new RegressionParams(min, max, detect, filter, rep, ifTemplate, corrLim, 1, regParams)
    }
  }
  def estimatePeaks(features: PeptideFeatures, regParams: RegressionParams): Int = {
    return (features.average * regParams.regressionParams(0) + features.lowpassaverage * regParams.regressionParams(1) + features.fakelowpassaverage * regParams.regressionParams(2) + features.max * regParams.regressionParams(3) + features.median * regParams.regressionParams(4) + features.numberOfChanges * regParams.regressionParams(5) + features.variance * regParams.regressionParams(6)).toInt

  }
}