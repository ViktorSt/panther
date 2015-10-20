package se.lth.immun
import akka.actor._
import akka.actor.{ Actor, ActorRef, Props }
import java.net.InetSocketAddress
import se.lth.immun.protocol.MSDataProtocolActors
import se.lth.immun.protocol.MSDataProtocol._



/**
 * @author viktor
 */
object WaitForPanther {
  
  var infinitePoller: ActorRef = _
  val tomte = GetStatus.newBuilder
  val request = MasterRequest.newBuilder.setGetStatus(tomte).build
  
  def main(args: Array[String]) = {
    println("Starting")
    
    val server = new InetSocketAddress("localhost", 12345)
    val system = ActorSystem()
    
    infinitePoller = system.actorOf(Props[TramlPoller])
    val client = system.actorOf(MSDataProtocolActors.ClientInitiator.props(server, infinitePoller))
    system.awaitTermination()
  }
  class TramlPoller extends Actor {

    import MSDataProtocolActors._

    def receive = {
      case msg: String =>
        println(msg)

      case MSDataProtocolConnected(remote, local) =>
        println("connected")
          sender ! request

      case MSDataReply(msg, nBytes, checkSum, timeTaken, remote) =>
        println("data recieved")
        println(msg.getStatus.toString())
        if(msg.getStatus.getProgress==msg.getStatus.getProgressMax){
          println("closing")
          System.exit(0)
        }else{
          Thread.sleep(10000)
          sender ! request
        }

    }
  }
}