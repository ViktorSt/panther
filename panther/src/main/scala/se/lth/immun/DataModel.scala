package se.lth.immun

/**
 * @author viktor
 */
object DataModel {
  
  val SMALL_ENOUGH_RT_DIFF = 1.0
  
  case class FragmentArea(fragment:String, area:Double)
  
  class Peak(
      val sequence:String,
      val iStart: Int,
      val iApex: Int,
      val iEnd: Int,
      val rtStart:Double,
      val rtApex:Double,
      val rtEnd:Double,
      val fragmentAreas:Seq[FragmentArea],
      var truePeak: Boolean
  ) {
    
    def isTheSame(p:Peak):Boolean = {
      math.abs(rtApex - p.rtApex) < SMALL_ENOUGH_RT_DIFF &&
        fragmentAreas.exists(fa1 => 
          p.fragmentAreas.exists(fa2 => 
            fa1.fragment == fa2.fragment && closeEnough(fa1.area, fa2.area)
           )
          )
    }
    
    def closeEnough(area1:Double, area2:Double) =
      area1 > area2 * 0.1 && area2 > area1 * 0.1
  }
}