#include "ClientSocket.h"
#include "SocketException.h"
#include <iostream>
#include <string>

int main ( int argc, char *argv[] )
{
  try
    {

      ClientSocket client_socket ( "localhost", 32000 );

      std::string reply;
      NN_to_tennix msg;
      tennix_to_NN resp;

      //d.i[0] = 33;
      //d.i[1] = 44;

      try
	{
	  client_socket << msg;//"Test message.";
	  client_socket >> resp;
          std::cout << "Response seq no :" << resp.msg.seq_no << " | " <<resp.msg.ball_x << std::endl;

	}
      catch ( SocketException& ) {}

      //std::cout << "We received this response from the server:\n\"" << d.i <<"\"\n";;

    }
  catch ( SocketException& e )
    {
      std::cout << "Exception was caught:" << e.description() << "\n";
    }

  return 0;
}
