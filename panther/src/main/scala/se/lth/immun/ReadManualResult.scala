package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.BufferedReader
import java.io.BufferedWriter
import java.io.BufferedReader

import java.io.FileWriter;
import java.util.Scanner
/**
 * @author viktor
 */
object ReadManualResult {
  def main(args: Array[String]) = {
    var name = "napedro_L120227_008_SW"

    val filename = new File(args(0))
    if (args.length == 2) {

      name = args(1)
    }
    val output = new File("result")
    val bufferedReader = new BufferedReader(new FileReader(filename))
    val bufferedWriter = new BufferedWriter(new FileWriter(output))
    var temp = bufferedReader.readLine()
    bufferedWriter.write(name)
    bufferedWriter.newLine()
    while (temp != null) {
      val temp2 = temp.split(';')
      if (temp2(3) == name) {
        bufferedWriter.write(temp2(temp2.size - 1) + "\t" + temp2(1) + "\t " + temp2(8) + "\t" + temp2(9))
        bufferedWriter.newLine()
      }
      temp = bufferedReader.readLine()
    }
    bufferedWriter.close()
  }
}