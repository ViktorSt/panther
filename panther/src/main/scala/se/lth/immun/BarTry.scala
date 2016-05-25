package se.lth.immun
import scalax.chart._
import scalax.chart.module.ChartFactories
import java.io.BufferedReader
import java.io.FileReader
/**
 * @author viktor
 */
object BarTry {
  def main(args: Array[String]) = {
    val bufferedReader = new BufferedReader(new FileReader("data"))
    var temp = bufferedReader.readLine()
    var findingRate = 0.0
    var falsePositiveRate = 0.0
    var count = 0
    var filename = ""
    var min = 0
    var max = 0
    var peakDetectionAndGrouping = 0
    var filter = 0
    var time = 0.0
    var i = 0
    val ds = new org.jfree.data.category.DefaultCategoryDataset
    while (temp != null && i < 100) {
      count = count % 9
      count match {
        case 0 => filename = temp
        case 1 =>
          min = temp.split(": ")(1).split(" boundaryMax")(0).toInt
          max = temp.split(": ")(2).toInt
        case 2 => peakDetectionAndGrouping = temp.split(": ")(1).toInt
        case 3 => filter = temp.split(": ")(1).toInt
        case 4 =>
        case 5 => findingRate = temp.split(": ")(1).toDouble
        case 6 => falsePositiveRate = temp.split(": ")(1).toDouble
        case 7 => time = temp.split(": ")(1).split("ms")(0).toDouble
        case 8 =>
        ds.addValue(findingRate, min, max)
//        ds.addValue(falsePositiveRate, peakDetectionAndGrouping, "falsePositiveRate")
      }
      temp = bufferedReader.readLine()
      i+=1
      
      count += 1
    }
    val chart = ChartFactories.BarChart(ds)
    chart.show()

  }
}