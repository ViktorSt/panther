package se.lth.immun

/**
 * @author viktor
 */
class FragFeatures(y: Seq[Double], name: String) {
  var fragName=name
  var average = 0.0
  var variance = 0.0
  var median=0.0
  var max = 0.0
  var numberOfChanges=0
  var lowpassaverage=0.0
  var fakelowpassaverage=0.0
  var counter =0
    val dy=getDerivate(y)
    var averagePower2 =0.0
    for(i<-0 until y.size){
      if(i+1<dy.size)
      if(math.signum(dy(i)) != math.signum(dy(i + 1))){
        numberOfChanges+=1
      }
      average+=y(i)/y.size
      averagePower2+=math.pow(y(i),2)/y.size
      if(y(i)>max){
        max=y(i)
      }
      if(y(i)>50){
        lowpassaverage+=y(i)
        counter+=1
      }
    }
  if(counter!=0){
  fakelowpassaverage=lowpassaverage/counter
  }
  lowpassaverage=lowpassaverage/y.size    
  variance= averagePower2-math.pow(average,2)
      
    val a =y.sorted
    median=(a(a.size/2)+a((a.size+1)/2))/2
  
  def getDerivate(x: Seq[Double]): Seq[Double] =
    return x.zip(x.tail).map(tu => tu._2 - tu._1)
}