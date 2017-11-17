#ifndef TOKENBUCKET_HH
#define TOKENBUCKET_HH

#include <vector>
#include <deque>

class Token_Bucket
{
private:
    double _token_rate;
    std::vector< double > _token;
    std::vector< double > _last_tickno;
    
    double _max_token;

public:
    Token_Bucket(const double token_rate = 10,double max_token = 100):_token_rate(token_rate),_token(),_last_tickno(),_max_token(max_token){};

    void autosize(const unsigned int id){
        if(_token.size() <= id){
            _token.resize( id + 1 );
            _last_tickno.resize( id + 1 );
        }
    }
    void set_rate(const double rate){
        _token_rate = rate;
    }
    void set_max_token(const unsigned int max_token){
        _max_token = max_token;
    }
    double get_token(const unsigned int src ,const double &tickno){
        autosize(src);

        double delta = tickno - _last_tickno[src];
        _token[src] += delta*_token_rate;
        _token[src] = _token[src]>_max_token?_max_token:_token[src];
        _last_tickno[src] = tickno;
        return _token[src];
    }
    bool used_token(const unsigned int src,const double &tickno){
        double tk = get_token(src,tickno);
        if(tk < 1) return false;
        _token[src] -= 1;
        return true;
    }
    void reset(double token_rate = -1)
    {
        _token.clear();
        _last_tickno.clear();
        if(token_rate > 0)
            _token_rate = token_rate;
    }
};

#endif
