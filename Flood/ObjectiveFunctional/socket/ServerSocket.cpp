// Implementation of the ServerSocket class

#include "ServerSocket.h"
#include "SocketException.h"


ServerSocket::ServerSocket ( int port )
{
  if ( ! Socket::create() )
    {
      throw SocketException ( "Could not create server socket." );
    }

  if ( ! Socket::bind ( port ) )
    {
      throw SocketException ( "Could not bind to port." );
    }

  if ( ! Socket::listen() )
    {
      throw SocketException ( "Could not listen to socket." );
    }

}

ServerSocket::~ServerSocket()
{
}


const ServerSocket& ServerSocket::operator << ( const std::string& s ) const
{
  if ( ! Socket::send ( s ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;

}



const ServerSocket& ServerSocket::operator >> ( std::string& s ) const
{
  if ( ! Socket::recv ( s ) )
    {
      throw SocketException ( "Could not read from socket." );
    }

  return *this;
}

/*
 * Test to send struct
 */

const ServerSocket& ServerSocket::operator << ( tennix_to_NN& d ) const
{
  if ( ! Socket::send ( d ) )
    {
      throw SocketException ( "Could not write to socket." );
    }

  return *this;
}

const ServerSocket& ServerSocket::operator >> ( NN_to_tennix& d ) const
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

void ServerSocket::accept ( ServerSocket& sock )
{
  if ( ! Socket::accept ( sock ) )
    {
      throw SocketException ( "Could not accept socket." );
    }
}
