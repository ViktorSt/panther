package se.lth.immun


/**
 * @author viktor
 */
class RegressionParams(minimum: Int, maximum: Int, detection: Int, filt: Int, repetitions: Int, ifTemplate: Int, correlationLimit: Double, thresholdOn:Int, regParams: Array[Double]) {
  val min = minimum
  val max = maximum
  val detect = detection
  val filter = filt
  val rep = repetitions
  val ifTemp = ifTemplate
  val corrLim = correlationLimit
  val ifThresholdOn= thresholdOn
  val regressionParams = regParams
 
}
