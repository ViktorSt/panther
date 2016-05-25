package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import se.lth.immun.FromTramToTable.Trace

/**
 * @author viktor
 */
class PeptideFeatures(name: String) {
  var peptideName = name
  var average = 0.0
  var variance = 0.0
  var median = 0.0
  var max = 0.0
  var numberOfChanges = 0
  var peakFound = false
  var fragments = new ArrayBuffer[FragFeatures]
  var fakelowpassaverage=0.0
  var lowpassaverage=0.0

 
}