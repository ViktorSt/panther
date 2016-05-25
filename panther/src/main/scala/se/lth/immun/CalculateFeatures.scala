package se.lth.immun

import scala.collection.mutable.ArrayBuffer

import se.lth.immun.FromTramToTable.Trace
/**
 * @author viktor
 */
class CalculateFeatures {
  def calcFeatures(database: ArrayBuffer[Trace], traceIds: ArrayBuffer[String]): ArrayBuffer[PeptideFeatures] = {
    var rt = new ArrayBuffer[PeptideFeatures]
    var old = ""
    var j = -1
    for (i <- 0 until database.size) {
      var y = database(i).intensities
      val ids = traceIds(i).split('/')
      val peptideId = ids(0) + "/" + ids(1).charAt(0)
      val fragId = ids(1).split('_')(1)
      if (old == peptideId) {
        rt(j).fragments += new FragFeatures(y, fragId)
      }
      if (peptideId != old || database.size - 1 == i) {
        if (j > -1) {
          var tempMax = 0.0
          val medians = rt(j).fragments.map(x => x.median).sorted
          rt(j).median = (medians(medians.size / 2) + medians((medians.size + 1) / 2)) / 2
          for (k <- 0 until rt(j).fragments.size) {
            if (rt(j).fragments(k).max > tempMax) {
              tempMax = rt(j).fragments(k).max

            }
            rt(j).average += rt(j).fragments(k).average / rt(j).fragments.size
            rt(j).lowpassaverage += rt(j).fragments(k).lowpassaverage / rt(j).fragments.size
            rt(j).fakelowpassaverage+=rt(j).fragments(k).fakelowpassaverage/rt(j).fragments.size
            rt(j).numberOfChanges += rt(j).fragments(k).numberOfChanges
            rt(j).variance += rt(j).fragments(k).variance

          }
          rt(j).max = tempMax
        }
        if (i != database.size - 1) {
          j += 1
          rt += new PeptideFeatures(peptideId)
          rt(j).fragments += new FragFeatures(y, fragId)
          old = peptideId
        }
      }

    }
    return rt
  }
}