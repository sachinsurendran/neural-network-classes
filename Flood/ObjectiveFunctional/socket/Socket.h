// Definition of the Socket class

#ifndef Socket_class
#define Socket_class


#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <string>
#include <arpa/inet.h>
#include <iostream>

/*
 *  Test struct send
 */

struct NN_to_tennix_msg {
	int msg_type; /* Type of message */
	int keys[3];  /* UP, DOWN, HIT */
	int seq_no;
};

struct tennix_to_NN_msg {
	int msg_type;
	float opponent_x;
	float opponent_y;
	float ball_x;
	float ball_y;
	int seq_no;
        float fitness;
};


class NN_to_tennix {
    public:
        //int i[10];
        struct NN_to_tennix_msg msg;
        NN_to_tennix(void);
        int size(void);
};

class tennix_to_NN {
    public:
        struct tennix_to_NN_msg msg;
        tennix_to_NN(void);
        int size(void);
};


const int MAXHOSTNAME = 200;
const int MAXCONNECTIONS = 5;
const int MAXRECV = 500;

class Socket
{
 public:
  Socket();
  virtual ~Socket();

  // Server initialization
  bool create();
  bool bind ( const int port );
  bool listen() const;
  bool accept ( Socket& ) const;

  // Client initialization
  bool connect ( const std::string host, const int port );

  // Data Transimission
  bool send ( const std::string ) const;
  int recv ( std::string& ) const;

  bool send ( const NN_to_tennix ) const;
  int recv ( tennix_to_NN& ) const;
  bool send ( const tennix_to_NN) const;
  int recv ( NN_to_tennix& ) const;


  void set_non_blocking ( const bool );

  bool is_valid() const { return m_sock != -1; }

 private:

  int m_sock;
  sockaddr_in m_addr;


};




#endif
