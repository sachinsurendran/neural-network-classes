// Definition of the ClientSocket class

#ifndef ClientSocket_class
#define ClientSocket_class

#include "Socket.h"


class ClientSocket : private Socket
{
 public:

  //ClientSocket ( std::string host, int port );
  //Introduced to get around the initialization when client socket is declared in another class
  void Initialize (std::string host, int port );
  virtual ~ClientSocket(){};

  const ClientSocket& operator << ( const std::string& ) const;
  const ClientSocket& operator >> ( std::string& ) const;

  const ClientSocket& operator << ( NN_to_tennix& d ) const;
  const ClientSocket& operator >> ( tennix_to_NN& d ) const;

};


#endif
