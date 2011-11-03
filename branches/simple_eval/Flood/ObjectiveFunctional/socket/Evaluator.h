#ifndef __EVALUATOR_H__
#define __EVALUATOR_H__

#include "ClientSocket.h"
#include "SocketException.h"
#include <iostream>
#include <string>
#include <cstring>

#define KEY_LEN 3 * sizeof(int)

#define TENNIX_SERVER_PORT 32000

#define TRUE  1
#define FALSE 0
#define ERROR -1


enum key_types {
	UP   = 0,
	DOWN = 1,
	HIT  = 2
};

enum msg_type {
    GAME_INIT = 0,
    NN_RESPONSE = 1,
    END_GAME  = 2
};

class GameState
{
    public:
        float darwin_x;
        float darwin_y;
	float opponent_x;
	float opponent_y;
	float ball_x;
	float ball_y;
        float fitness;
        char  game_ended;/* Bool which tell game ending, time to check fitness */
};

class Evaluator
{
    private:

        ClientSocket client_socket;
        int seq_no;
        int initiate_evaluation_session();

    public:
        Evaluator();
        int get_game_state(GameState *gamestate);
        int send_NN_response(int *keys, int len);
};

int
Evaluator::initiate_evaluation_session(void)
{
    NN_to_tennix init_msg;

    init_msg.msg.msg_type = GAME_INIT;

      try
      {
	client_socket << init_msg;
      }
      catch ( SocketException& ) 
      {
          std::cout << "Error in" << __FUNCTION__ << std::endl;
          return ERROR;
      }
}


Evaluator::Evaluator() {

    seq_no = 0;

    try
    {
      client_socket.Initialize( "localhost", 32000 );
    }
    catch ( SocketException& e )
    {
      std::cout << "Exception was caught:" << e.description() << "\n";
    }

    initiate_evaluation_session();
}

int
Evaluator::get_game_state(GameState *gamestate)
{
    tennix_to_NN tennix_resp_msg;

    /* Get the response from tennix */
    client_socket >> tennix_resp_msg;

    //std::cout << "Darwin  X = " << tennix_resp_msg.msg.darwin_x << "Y = " << tennix_resp_msg.msg.darwin_y << std::endl; 

    gamestate->darwin_x = tennix_resp_msg.msg.darwin_x;
    gamestate->darwin_y = tennix_resp_msg.msg.darwin_y;
    gamestate->opponent_y = tennix_resp_msg.msg.opponent_y;
    gamestate->opponent_x = tennix_resp_msg.msg.opponent_x;
    gamestate->ball_y = tennix_resp_msg.msg.ball_y;
    gamestate->ball_x = tennix_resp_msg.msg.ball_x;
    gamestate->fitness= tennix_resp_msg.msg.fitness;

    if (tennix_resp_msg.msg.msg_type == END_GAME)
    {
        gamestate->game_ended = true;
    }
}

int 
Evaluator::send_NN_response(int *keys, int len)
{
    NN_to_tennix NN_resp_msg;

    if ( len != (KEY_LEN / sizeof(int)))
        std::cout << "Error: Key length incorrect "<< std::endl;

    NN_resp_msg.msg.msg_type = NN_RESPONSE;
    memcpy(NN_resp_msg.msg.keys, keys, KEY_LEN);
    NN_resp_msg.msg.seq_no = seq_no++;

    /* Send the response */
    client_socket << NN_resp_msg;
}

#endif
