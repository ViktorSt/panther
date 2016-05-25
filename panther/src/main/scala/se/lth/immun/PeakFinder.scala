package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import se.lth.immun.FromTramToTable.PC
import se.lth.immun.FromTramToTable.RTPc
import se.lth.immun.FromTramToTable.PeakAndFrags

/**
 * @author viktor
 */
object PeakFinder {
  trait SinglePeakFinder{
    def findPc(y: Seq[Double], rts: Seq[Double]): Seq[PC]    
  }
  
  case class SingleDerivativeDetect(boundaryMin: Int, boundaryMax:Int) extends SinglePeakFinder{
    def  findPc(y: Seq[Double], rts: Seq[Double]): Seq[PC]  =
      bruteForce(y, rts, boundaryMin, boundaryMax)
  }
  case class TemplateDetect(correlationLimit: Double) extends SinglePeakFinder{
        def  findPc(y: Seq[Double], rts: Seq[Double]): Seq[PC]  =
          findTemplatePCs(y, rts, correlationLimit)
  }
  case class MultiDerivativeDetect(boundaryMin: Int,boundaryMax: Int){
    def findPc(ys: ArrayBuffer[Seq[Double]], listOfFrags: ArrayBuffer[String], rts: Seq[Double]): Seq[PeakAndFrags] = 
      findMultiPeaks(ys, listOfFrags, rts, boundaryMin, boundaryMax)
  }
  val areaFactor=1.771
  def getDerivate(x: Seq[Double]): Seq[Double] =
    return x.zip(x.tail).map(tu => tu._2 - tu._1)
    
  def bruteForce(
    y: Seq[Double],
    rts: Seq[Double],
   
    boundaryMin: Int,
     boundaryMax: Int): Seq[PC] = {
    var ret = new ArrayBuffer[PC]
    
    val dy =  getDerivate(y)
    val ddy = getDerivate(dy)
    var istart = 0
    var iapex = -1
    for (i <- 0 until dy.length - 1) {
      val starty = y(istart)
      if (math.signum(dy(i)) != math.signum(dy(i + 1))) {
        if (ddy(i) > 0) { // local min
          if (iapex >= 0 && boundaryMin >= y(i + 1) - starty) {
            //            println("valid min")
            var max = y(iapex)
            val endy = y(i + 1)
            if (boundaryMax <= max - math.max(starty, endy)) {
              //              println("found peak")
              var area = 0.0
              for (j <- istart until i + 1) {
                if (j == 0) {

                  area += y(j) * (rts(j)) * 60
                } else {
                  area += y(j) * (rts(j) - rts(j - 1)) * 60
                }
              }
              ret += new PC(istart, iapex, i + 1, new ArrayBuffer[RTPc] += new RTPc(rts(istart), rts(iapex), rts(i + 1), area* areaFactor))
            }

            istart = i + 1

            iapex = -1
          } else if (y(i + 1) < starty) {
            istart = i + 1
          }
        } else { // local max
          if (iapex == (-1)) {
            //            println("setting max, old apex -1")
            iapex = i + 1
          } else if (y(iapex) < y(i + 1)) {

            //            println("setting max, new apex higher")
            iapex = i + 1
          }

        }
      } else if (y(i) <= starty && i == istart + 1) { //concider if best approach
        istart = i
      }
    }

    return ret
  }
  
  val template = Array(125.1, 469.8, 1273.8, 2340.6, 3149.2, 2852.1, 1694.4, 685.6, 280.3);
   def findTemplatePCs(y: Seq[Double], rts: Seq[Double], correlationLimit: Double): Seq[PC] = {
    var ret = new ArrayBuffer[PC]
    for (i <- 0 until y.length - template.length) {

      var xiyi = 0.0
      var xi = 0.0
      var yi = 0.0
      var xi2 = 0.0
      var yi2 = 0.0
      var area = 0.0
      var iapex = -1
      for (j <- 0 until template.length) {
        if (iapex == -1) {
          iapex = i + j
        } else if (y(i + j) > y(iapex)) {
          iapex = i + j
        }
        xiyi += template(j) * y(i + j)
        xi += y(i + j)
        yi += template(j)
        xi2 += math.pow(y(i + j), 2)
        yi2 += math.pow(template(j), 2)
        if (i + j == 0) {

          area += y(i + j) * (rts(i + j)) * 60
        } else {
          area += y(i + j) * (rts(i + j) - rts(i + j - 1)) * 60
        }
      }
      val rxy = (9 * xiyi - xi * yi) / (math.sqrt(9 * xi2 - math.pow(xi, 2)) * math.sqrt(9 * yi2 - math.pow(yi, 2)))
      if (rxy > correlationLimit) {

        ret += new PC(i, iapex, i + template.length - 1, new ArrayBuffer[RTPc] += new RTPc(rts(i), rts(iapex), rts(i + template.length - 1), area *areaFactor))
        //println(rts(i) + " " + rts(i + template.length - 1) + " " + rxy)
      }
    }

    return ret
  }
   def findMultiPeaks(ys: ArrayBuffer[Seq[Double]], listOfFrags: ArrayBuffer[String], rts: Seq[Double], boundaryMin: Int, boundaryMax: Int): Seq[PeakAndFrags] = {
    var dys = new ArrayBuffer[Seq[Double]]
    var ddys = new ArrayBuffer[Seq[Double]]
    var ret = new ArrayBuffer[PeakAndFrags]
    var minflag = false
    var maxflag = false

    var flags = Array(false, false, false, false)
    for (pl <- ys) {
      val tempDeriv = getDerivate(pl)
      dys += tempDeriv
      ddys += getDerivate(tempDeriv)
    }
    var istart = 0
    var iapex = -1
    var starty = ys(0)(istart)
    var endy = -1.0
    var max = -1.0
    for (i <- 0 until dys(0).length - 1) {
      var temprtPc = new ArrayBuffer[RTPc]
      var j = -1
      var k = 0
      var z = -1
      var l = 0
      for (t <- 0 until dys.length) {
        j += 1
        if (i == 0) {
          if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
            if (ddys(t)(i + 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        } else if (i == dys(0).length - 2) {
          if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
            if (ddys(t)(i - 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        } else {
          if (math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
            if (ddys(t)(i + 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
            if (ddys(t)(i - 1) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          } else if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1))) {
            if (ddys(t)(i) > 0) {

              k += 1
              z = j
            } else {
              l += 1
            }
          }
        }

      }

      if (l > 1) { // local max
        var tempMax = ys(0)(i + 1)
        for (y <- ys) {
          if (y(i + 1) > tempMax) {
            tempMax = y(i + 1)
          }
        }
        if (iapex == -1) {
          iapex = i + 1
          max = tempMax
        } else if (tempMax > max) {
          iapex = i + 1
          max = tempMax
        }
        for (t <- 0 until dys.length) {
          if (i == 0) {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2))) {
              flags(t) = true
            }
          } else if (i == dys(0).length - 2) {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
              flags(t) = true
            }
          } else {
            if (math.signum(dys(t)(i)) != math.signum(dys(t)(i + 1)) || math.signum(dys(t)(i + 1)) != math.signum(dys(t)(i + 2)) || math.signum(dys(t)(i - 1)) != math.signum(dys(t)(i))) {
              flags(t) = true
            }
          }
        }

      } else if (k > 1) { // local min

        for (y <- ys) {
          if (boundaryMin >= y(i + 1) - starty) {
            minflag = true
          }
        }
        if (iapex >= 0 && minflag) {
          var tempMin = ys(0)(i + 1)
          for (y <- ys) {
            if (y(i + 1) < tempMin) {
              tempMin = y(i + 1)
            }
          }
          endy = tempMin
          if (boundaryMax <= max - math.max(starty, endy)) {
            var tempFrags = new ArrayBuffer[String]
            for (t <- 0 until ddys.length) {
              if (flags(t)) {
                tempFrags += listOfFrags(t)
                var tempArea = 0.0
                for (j <- istart until i + 1) {
                  if (j == 0) {

                    tempArea += ys(t)(j) * (rts(j)) * 60
                  } else {
                    tempArea += ys(t)(j) * (rts(j) - rts(j - 1)) * 60
                  }
                }
                temprtPc += new RTPc(rts(istart), rts(iapex), rts(i + 1), tempArea * areaFactor)
              }
            }
            ret += new PeakAndFrags(new PC(istart, iapex, i + 1, temprtPc), tempFrags)
          }

          flags = Array(false, false, false, false)
          iapex = -1

          istart = i + 1
        }

        minflag = false

      } else if (ys(0)(i) < ys(0)(istart) && i == istart + 1) {
        istart = i
      } else {
      }
    }
    return ret
  }
}