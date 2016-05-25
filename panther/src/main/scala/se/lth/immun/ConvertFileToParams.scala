package se.lth.immun

import java.util.Scanner
import java.io.File

/**
 * @author viktor
 */
object ConvertFileToParams {
  def main(args: Array[String]) = {
    println("STARTING")
    val fileReader = new Scanner(new File("coef.txt"))
    
       val fileWriter = new java.io.PrintWriter(new java.io.FileWriter("coefFile", true))
    while (fileReader.hasNext()) {
      println("got line")
      val line = fileReader.nextLine()
      val splitData = line.split(" ")
      val tempData = splitData(0)
      val usefulData = tempData.slice(7, tempData.length)
      val paramInt = Integer.parseInt(usefulData)
      
        var tempString=""
      if(paramInt<865){
        
        val dg=(paramInt-1)/288
        
        var temp = (paramInt-1)%288
        
        val filter=temp/36
        
        temp=temp%36
        val min = temp/12
        temp=temp%12
        
        val max = temp/4
        temp=temp%4
        val rep=temp
       min match {
          case 0 => tempString+="0(0) "
          case 1 => tempString+="30(1) "
          case 2 => tempString+="60(2) "
        }
       max match {
          case 0 => tempString+="0(0) "
          case 1 => tempString+="30(1) "
          case 2 => tempString+="60(2) "
        }
       dg match {
          case 0 => tempString+="BruteForce(1) "
          case 1 => tempString+="Sorted(2) "
          case 2 => tempString+="Multi(3) "
        }
       filter match {
          case 0 => tempString+="NoFilter(0) "
          case 1 => tempString+="SavitzkyGolay(1) "
          case 2 => tempString+="Daubechies2(2) "
          case 3 => tempString+="Daubechies3(3) "
          case 4 => tempString+="Daubechies8(4) "
          case 5 => tempString+="Haar(5) "
          case 6 => tempString+="ChainWavelet(6) "
          case 7 => tempString+="Multiwavelet(7) "
        }
       tempString+= rep+"("+rep+") NoTemplate(0) N/A(0)  "
       fileWriter.print(tempString)
       for(i<-1 until splitData.length){
         fileWriter.print(splitData(i)+ " ")
       }
       fileWriter.println()
      }
      else{
        var temp= paramInt-865
        
        val dg=temp/128
        
        temp =temp%128
        
        val filter=temp/16
        temp=temp%16
        
        val rep=temp/4
        temp =temp%4
        val korrlim= temp%4
        tempString+="N/A(0) N/A(0) "
        dg match {
          case 0 => tempString+="BruteForce(1) "
          case 1 => tempString+="Sorted(2) "
        }
        filter match {
          case 0 => tempString+="NoFilter(0) "
          case 1 => tempString+="SavitzkyGolay(1) "
          case 2 => tempString+="Daubechies2(2) "
          case 3 => tempString+="Daubechies3(3) "
          case 4 => tempString+="Daubechies8(4) "
          case 5 => tempString+="Haar(5) "
          case 6 => tempString+="ChainWavelet(6) "
          case 7 => tempString+="Multiwavelet(7) "
        }
        tempString +=rep + "(" +rep + ") Template(1) "
        korrlim match {
          case 0 => tempString+="0.49(0.49) "
          case 1 => tempString+="0.74(0.74) "
          case 2 => tempString+="0.95(0.95) "
          case 3 => tempString+="0.99(0.99)"
        }
       fileWriter.print(tempString)
       for(i<-1 until splitData.length){
         fileWriter.print(splitData(i)+ " ")
       }
       fileWriter.println()
      }
    }
    fileWriter.flush()
    println("done")
  }
}