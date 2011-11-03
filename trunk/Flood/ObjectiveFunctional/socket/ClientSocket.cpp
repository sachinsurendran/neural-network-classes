// Implementation of the ClientSocket class

#include "ClientSocket.h"
#include "SocketException.h"
/*
ClientSocket::ClientSocket ( std::string host, int port )
{
  if ( ! Socket::create() )
    {
      throw SocketException ( "Could not create client socket." );
    }

  if ( ! Socket::connect ( host, port ) )
    {
      throw SocketException ( "Could not bind to port." );
    }

}
*/
void ClientSocket::Initialize ( std::string host, int port )
{
  if ( ! Socket::create() )
    {
      throw SocketException ( "Could not create client socket." );
    }

  /* Commenting this part as the Tennix server does not have code
   * to accept connection */

  if ( ! Socket::connect ( host, port ) )
    {
      throw SocketException ( "Could not bind to port." );
    }

}


const ClientSocket& ClientSocket::operator << ( const std::string& s ) const
{
  if ( ! Socket::send ( s ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;

}


const ClientSocket& ClientSocket::operator >> ( std::string& s ) const
{
  if ( ! Socket::recv ( s ) )
    {
      throw SocketException ( "Could not read from socket." );
    }

  return *this;
}

/*
 * Test to verify structure send
 */

const ClientSocket& ClientSocket::operator << ( NN_to_tennix& d ) const
{
  if ( ! Socket::send ( d ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;

}


const ClientSocket& ClientSocket::operator >> ( tennix_to_NN& d ) const
{
  if ( ! Socket::recv ( d ) )
    {
      throw SocketException ( "Could not read from socket." );
    }

  return *this;
}

/*
 * End of test
 */
