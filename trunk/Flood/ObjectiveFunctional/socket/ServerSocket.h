// Definition of the ServerSocket class

#ifndef ServerSocket_class
#define ServerSocket_class

#include "Socket.h"


class ServerSocket : private Socket
{
 public:

  ServerSocket ( int port );
  ServerSocket (){};
  virtual ~ServerSocket();

  const ServerSocket& operator << ( const std::string& ) const;
  const ServerSocket& operator >> ( std::string& ) const;

  const ServerSocket& operator << ( tennix_to_NN& d ) const;
  const ServerSocket& operator >> ( NN_to_tennix& d ) const;

  void accept ( ServerSocket& );

};


#endif
