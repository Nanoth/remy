#include <cassert>
#include <limits>

#include "receiver.hh"

Receiver::Receiver()
  : _collector(),_bucket(),_bucket_flag(0)
{
}

Receiver::Receiver(const double token_rate):_collector(),_bucket(token_rate),_bucket_flag(1)
{
}

void Receiver::accept( const Packet & p, const double & tickno ) noexcept
{
  autosize( p.src );
    
  if(_bucket_flag && _bucket.used_token(p.src,tickno)){
     _collector[ p.src ].push_back( p );
     _collector[ p.src ].back().tick_received = tickno;
  }
}

void Receiver::autosize( const unsigned int index )
{
  if ( index >= _collector.size() ) {
    _collector.resize( index + 1 );
  }
}

double Receiver::next_event_time( const double & tickno ) const
{
  for ( const auto & x : _collector ) {
    if ( not x.empty() ) {
      return tickno;
    }
  }
  return std::numeric_limits<double>::max();
}
