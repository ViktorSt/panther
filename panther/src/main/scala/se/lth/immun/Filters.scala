package se.lth.immun

/**
 * @author viktor
 */

object Filters {
  
  trait Filter {
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double]
  }
  case object NoFilter extends Filter {
       def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      signal
  }
  case object SavitzkyGolay9 extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      savitzkyGolay9(signal)
    
  }
  case object Daubchies2 extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      daubchies2(signal, nReps, thresholdOn)
    
  }
  case object Daubchies3 extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      daubchies3(signal, nReps, thresholdOn)
    
  }
  case object Daubchies8 extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      daubchies8(signal, nReps, thresholdOn)
    
  }
  case object Haar extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      haar(signal, nReps, thresholdOn)
    
  }
  case object ChainWavelet extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      chainWavelet(signal, nReps, thresholdOn)
    
  }
  case object MultiWavelet extends Filter {
    
    def applyFilter(signal:Seq[Double], nReps:Int, thresholdOn: Boolean):Seq[Double] = 
      multiWavelet(signal, nReps, thresholdOn)
    
  }
  
  val haar = Array(
    0.7071067811865475,
    0.7071067811865475)

  val gaar = Array(
    0.7071067811865475,
    -0.7071067811865475)

  val G8 = Array(
    0.0018899503329007,
    0.0003029205145516,
    -0.0149522583367926,
    -0.0038087520140601,
    0.0491371796734768,
    0.0272190299168137,
    -0.0519458381078751,
    -0.3644418948359564,
    0.7771857516997478,
    -0.4813596512592012,
    -0.0612733590679088,
    0.1432942383510542,
    0.0076074873252848,
    -0.0316950878103452,
    -0.0005421323316355,
    0.0033824159513594)
  val H8 = Array(
    -0.0033824159513594,
    -0.0005421323316355,
    0.0316950878103452,
    0.0076074873252848,
    -0.1432942383510542,
    -0.0612733590679088,
    0.4813596512592012,
    0.7771857516997478,
    0.3644418948359564,
    -0.0519458381078751,
    -0.0272190299168137,
    0.0491371796734768,
    0.0038087520140601,
    -0.0149522583367926,
    -0.0003029205145516,
    0.0018899503329007)
  val H3 = Array(
    math.sqrt(2) * (1 + math.sqrt(10) + math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (5 + math.sqrt(10) + 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (10 - 2 * math.sqrt(10) + 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (10 - 2 * math.sqrt(10) - 2 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (5 + math.sqrt(10) - 3 * math.sqrt(5 + 2 * math.sqrt(10))) / 32,
    math.sqrt(2) * (1 + math.sqrt(10) - math.sqrt(5 + 2 * math.sqrt(10))) / 32)
  val G3 = Array(H3(5), -H3(4), H3(3), -H3(2), H3(1), -H3(0))

  val H2 = Array(
    (1 + math.sqrt(3)) / (4 * math.sqrt(2)),
    (3 + math.sqrt(3)) / (4 * math.sqrt(2)),
    (3 - math.sqrt(3)) / (4 * math.sqrt(2)),
    (1 - math.sqrt(3)) / (4 * math.sqrt(2)))

  val G2 = Array(H2(3), -H2(2), H2(1), -H2(0))
  def savitzkyGolay9(a: Seq[Double]): Seq[Double] = {

    if (a.length < 9)
      return a

    var retArr = new Array[Double](a.length)
    var j = a.length
    for (i <- 0 until 4) {
      j -= 1
      retArr(i) = a(i)
      retArr(j) = a(j)
    }

    for (i <- 4 until a.length - 4) {
      retArr(i) = (0.417 * a(i) +
        +0.315 * (a(i - 1) + a(i + 1))
        + 0.070 * (a(i - 2) + a(i + 2))
        - 0.128 * (a(i - 3) + a(i + 3))
        + 0.035 * (a(i - 4) + a(i + 4)))
    }
    return retArr map (d => math.max(0, d))
  }
  
  def chainWavelet(x: Seq[Double], levels: Int, thresholdOn: Boolean): Seq[Double] = {
    var yh1 = new Array[Double](x.size / 2)
    var yl1 = new Array[Double](x.size / 2)
    val threshhold = 0.1
    var temp1 = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform
    {
      yh1(i) = x((i * 2) % x.size) * gaar(1) + x((i * 2 + 1) % x.size) * gaar(0)
      if(!thresholdOn){
        yh1(i)=0
      }
      else if (yh1(i) < 10) {
        yh1(i) = 0
      }
      yl1(i) = x((i * 2) % x.size) * haar(1) + x((i * 2 + 1) % x.size) * haar(0)
    }

    var yh2 = new Array[Double](yl1.size / 2)
    var yl2 = new Array[Double](yl1.size / 2)

    var temp2 = new Array[Double](yl1.size)
    for (i <- 1 until yl1.size / 2) //forward Transform
    {
      yh2(i) = yl1((i * 2) % yl1.size) * G2(3) + yl1((i * 2 + 1) % yl1.size) * G2(2) + yl1((i * 2 + 2) % yl1.size) * G2(1) + yl1((i * 2 + 3) % yl1.size) * G2(0)
      if(!thresholdOn){
        yh2(i)=0
      }
      else if (yh2(i) < 10) {
        yh2(i) = 0
      }
      yl2(i) = yl1((i * 2) % yl1.size) * H2(3) + yl1((i * 2 + 1) % yl1.size) * H2(2) + yl1((i * 2 + 2) % yl1.size) * H2(1) + yl1((i * 2 + 3) % yl1.size) * H2(0)
    }

    if (levels > 0) {
      chainWavelet(yl2, levels - 1, thresholdOn)
    }
    for (i <- 0 until yl1.size / 2) //inverse Transform
    {
      temp2((i * 2) % x.size) = yl2(((i - 1) % yl2.size + yl2.size) % yl2.size) * H2(1) + yl2((i) % yl2.size) * H2(3) + yh2((((yh2.size - 1 + i) % yh2.size) + yh2.size) % yh2.size) * G2(1) + yh2((yh2.size + i) % yh2.size) * G2(3)
      temp2((i * 2 + 1) % x.size) = yl2(((yl2.size - 1 + i) % yl2.size + yl2.size) % yl2.size) * H2(0) + yl2((yl2.size + i) % yl2.size) * H2(2) + yh2(((yh2.size - 1 + i) % yh2.size + yh2.size) % yh2.size) * G2(0) + yh2((yh2.size + i) % yh2.size) * G2(2)

    }

    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp1((i * 2) % x.size) = temp2(((i) % temp2.size + temp2.size) % temp2.size) * haar(1) + yh1((((yh1.size - 1 + i) % yh1.size) + yh1.size) % yh1.size) * gaar(1)
      temp1((i * 2 + 1) % x.size) = temp2(((temp2.size + i) % temp2.size + temp2.size) % temp2.size) * haar(0) + yh1(((yh1.size + i) % yh1.size + yh1.size) % yh1.size) * gaar(0)

    }

    return temp1.toSeq map (d => math.max(0, d))
  }

  def multiWavelet(x: Seq[Double], levels: Int, thresholdOn: Boolean): Seq[Double] = {
    var temp = daubchies2(x, 0, thresholdOn)
    temp = daubchies8(temp, 0, thresholdOn)
    if (levels > 0) {
      temp = multiWavelet(temp, levels - 1, thresholdOn)
    }

    return temp
  }

  def haar(x: Seq[Double], levels: Int, thresholdOn: Boolean): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    val threshhold = 0.1
    var temp = new Array[Double](x.size)
    yh(0) = 0
    yl(0) = x(0)
    for (i <- 1 until x.size / 2) //forward Transform why is i not 0?
    {
      yh(i) = x((i * 2) % x.size) * gaar(1) + x((i * 2 + 1) % x.size) * gaar(0)
      if(!thresholdOn){
        yh(i)=0
      }
      else if (yh(i) < 10) {
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * haar(1) + x((i * 2 + 1) % x.size) * haar(0)
    }
    if (levels > 0) {
      yl = haar(yl.toSeq, levels - 1 , thresholdOn).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i) % yl.size + yl.size) % yl.size) * haar(1) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * gaar(1)
      temp((i * 2 + 1) % x.size) = yl(((yl.size + i) % yl.size + yl.size) % yl.size) * haar(0) + yh(((yh.size + i) % yh.size + yh.size) % yh.size) * gaar(0)

    }

    return temp.toSeq map (d => math.max(0, d))
  }

  def daubchies2(x: Seq[Double], levels: Int, thresholdOn: Boolean): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    yh(0) = 0
    yl(0) = x(0)
    for (i <- 1 until x.size / 2) //forward Transform
    {
      yh(i) = x((i * 2) % x.size) * G2(3) + x((i * 2 + 1) % x.size) * G2(2) + x((i * 2 + 2) % x.size) * G2(1) + x((i * 2 + 3) % x.size) * G2(0)
      if(!thresholdOn){
        yh(i)=0
      }
      else if (yh(i) < 10) {
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * H2(3) + x((i * 2 + 1) % x.size) * H2(2) + x((i * 2 + 2) % x.size) * H2(1) + x((i * 2 + 3) % x.size) * H2(0)
    }
    if (levels > 0) {
      yl = daubchies2(yl.toSeq, levels - 1, thresholdOn).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 1) % yl.size + yl.size) % yl.size) * H2(1) + yl((i) % yl.size) * H2(3) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G2(1) + yh((yh.size + i) % yh.size) * G2(3)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H2(0) + yl((yl.size + i) % yl.size) * H2(2) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G2(0) + yh((yh.size + i) % yh.size) * G2(2)

    }

    return temp.toSeq map (d => math.max(0, d))
  }

  def daubchies3(x: Seq[Double], levels: Int, thresholdOn: Boolean): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until 2) {
      yh(i) = 0
      yl(i) = x(i)
    }
    for (i <- 2 until x.size / 2) //forward Transform again, why not i=0?
    {
      
      yh(i)= x((i*2)%x.size)*G3(5) + x((i*2+1)%x.size)*G3(4) + x((i*2+2)%x.size)*G3(3) + x((i*2+3)%x.size)*G3(2)+ x((i*2+4)%x.size)*G3(1)+ x((i*2+5)%x.size)*G3(0)
      if(!thresholdOn){
        yh(i)=0
      }
      else if(yh(i)<10){
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * H3(5) + x((i * 2 + 1) % x.size) * H3(4) + x((i * 2 + 2) % x.size) * H3(3) + x((i * 2 + 3) % x.size) * H3(2) + x((i * 2 + 4) % x.size) * H3(1) + x((i * 2 + 5) % x.size) * H3(0)
    }
    if (levels > 0) {
      yl = daubchies3(yl.toSeq, levels - 1, thresholdOn).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 2) % yl.size + yl.size) % yl.size) * H3(1) + yl(((i - 1) % yl.size + yl.size) % yl.size) * H3(3) + yl((i) % yl.size) * H3(5) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G3(1) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G3(3) + yh((yh.size + i) % yh.size) * G3(5)
      temp((i * 2 + 1) % x.size) = yl(((yl.size - 2 + i) % yl.size + yl.size) % yl.size) * H3(0) + yl(((yl.size - 1 + i) % yl.size + yl.size) % yl.size) * H3(2) + yl((yl.size + i) % yl.size) * H3(4) + yh(((yh.size - 2 + i) % yh.size + yh.size) % yh.size) * G3(0) + yh(((yh.size - 1 + i) % yh.size + yh.size) % yh.size) * G3(2) + yh((yh.size + i) % yh.size) * G3(4)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
  def daubchies8(x: Seq[Double], levels: Int,thresholdOn: Boolean): Seq[Double] = {
    var yh = new Array[Double](x.size / 2)
    var yl = new Array[Double](x.size / 2)
    var temp = new Array[Double](x.size)
    for (i <- 0 until x.size / 2) //forward Transform cause here i = 0....
    {

      yh(i) = x((i * 2) % x.size) * G8(15) + x((i * 2 + 1) % x.size) * G8(14) + x((i * 2 + 2) % x.size) * G8(13) + x((i * 2 + 3) % x.size) * G8(12) + x((i * 2 + 4) % x.size) * G8(11) + x((i * 2 + 5) % x.size) * G8(10) + x((i * 2 + 6) % x.size) * G8(9) + x((i * 2 + 7) % x.size) * G8(8) + x((i * 2 + 8) % x.size) * G8(7) + x((i * 2 + 9) % x.size) * G8(6) + x((i * 2 + 10) % x.size) * G8(5) + x((i * 2 + 11) % x.size) * G8(4) + x((i * 2 + 12) % x.size) * G8(3) + x((i * 2 + 13) % x.size) * G8(2) + x((i * 2 + 14) % x.size) * G8(1) + x((i * 2 + 15) % x.size) * G8(0)
      if(!thresholdOn){
        yh(i)=0
      }
      else if (yh(i) < 10) {
        yh(i) = 0
      }
      yl(i) = x((i * 2) % x.size) * H8(15) + x((i * 2 + 1) % x.size) * H8(14) + x((i * 2 + 2) % x.size) * H8(13) + x((i * 2 + 3) % x.size) * H8(12) + x((i * 2 + 4) % x.size) * H8(11) + x((i * 2 + 5) % x.size) * H8(10) + x((i * 2 + 6) % x.size) * H8(9) + x((i * 2 + 7) % x.size) * H8(8) + x((i * 2 + 8) % x.size) * H8(7) + x((i * 2 + 9) % x.size) * H8(6) + x((i * 2 + 10) % x.size) * H8(5) + x((i * 2 + 11) % x.size) * H8(4) + x((i * 2 + 12) % x.size) * H8(3) + x((i * 2 + 13) % x.size) * H8(2) + x((i * 2 + 14) % x.size) * H8(1) + x((i * 2 + 15) % x.size) * H8(0)
    }
    if (levels > 0) {
      yl = daubchies8(yl.toSeq, levels - 1, thresholdOn).toArray
    }
    for (i <- 0 until x.size / 2) //inverse Transform
    {
      temp((i * 2) % x.size) = yl(((i - 7) % yl.size + yl.size) % yl.size) * H8(1) + yl(((i - 6) % yl.size + yl.size) % yl.size) * H8(3) + yl(((i - 5) % yl.size + yl.size) % yl.size) * H8(5) + yl(((i - 4) % yl.size + yl.size) % yl.size) * H8(7) + yl(((i - 3) % yl.size + yl.size) % yl.size) * H8(9) + yl(((i - 2) % yl.size + yl.size) % yl.size) * H8(11) + yl(((i - 1) % yl.size + yl.size) % yl.size) * H8(13) + yl(((i) % yl.size + yl.size) % yl.size) * H8(15) + yh((((yh.size - 7 + i) % yh.size) + yh.size) % yh.size) * G8(1) + yh((((yh.size - 6 + i) % yh.size) + yh.size) % yh.size) * G8(3) + yh((((yh.size - 5 + i) % yh.size) + yh.size) % yh.size) * G8(5) + yh((((yh.size - 4 + i) % yh.size) + yh.size) % yh.size) * G8(7) + yh((((yh.size - 3 + i) % yh.size) + yh.size) % yh.size) * G8(9) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G8(11) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G8(13) + yh((((yh.size + i) % yh.size) + yh.size) % yh.size) * G8(15)
      temp((i * 2 + 1) % x.size) = yl(((yl.size + i - 7) % yl.size + yl.size) % yl.size) * H8(0) + yl(((yl.size + i - 6) % yl.size + yl.size) % yl.size) * H8(2) + yl(((yl.size + i - 5) % yl.size + yl.size) % yl.size) * H8(4) + yl(((yl.size + i - 4) % yl.size + yl.size) % yl.size) * H8(6) + yl(((yl.size + i - 3) % yl.size + yl.size) % yl.size) * H8(8) + yl(((yl.size + i - 2) % yl.size + yl.size) % yl.size) * H8(10) + yl(((yl.size + i - 1) % yl.size + yl.size) % yl.size) * H8(12) + yl(((yl.size + i) % yl.size + yl.size) % yl.size) * H8(14) + yh((((yh.size - 7 + i) % yh.size) + yh.size) % yh.size) * G8(1) + yh((((yh.size - 6 + i) % yh.size) + yh.size) % yh.size) * G8(3) + yh((((yh.size - 5 + i) % yh.size) + yh.size) % yh.size) * G8(5) + yh((((yh.size - 4 + i) % yh.size) + yh.size) % yh.size) * G8(7) + yh((((yh.size - 3 + i) % yh.size) + yh.size) % yh.size) * G8(9) + yh((((yh.size - 2 + i) % yh.size) + yh.size) % yh.size) * G8(11) + yh((((yh.size - 1 + i) % yh.size) + yh.size) % yh.size) * G8(13) + yh((((yh.size + i) % yh.size) + yh.size) % yh.size) * G8(15)

    }

    return temp.toSeq map (d => math.max(0, d))
  }
}