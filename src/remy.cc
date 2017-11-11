#include <cstdio>
#include <vector>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <queue>

#include "ratbreeder.hh"
#include "dna.pb.h"
#include "configrange.hh"
using namespace std;

typedef unsigned long long u64;
typedef long long s64 ;
int remycc_shift = 20;
double eps = 1e-8;

void print_range( const Range & range, const string & name )
{
  printf( "Optimizing for %s over [%f : %f : %f]\n", name.c_str(),
    range.low, range.incr, range.high );
}

struct M{
    u64 _send_ewma;
    u64 _rec_ewma;
    u64 _rtt_ratio;
    u64 _slow_rec_ewma;
    M():_send_ewma(0),
        _rec_ewma(0),
        _rtt_ratio(0),
        _slow_rec_ewma(0){};
    M(struct M& tmp):_send_ewma(tmp._send_ewma),
                    _rec_ewma(tmp._rec_ewma),
                    _rtt_ratio(tmp._rtt_ratio),
                    _slow_rec_ewma(tmp._slow_rec_ewma){};
    M& operator=(const Memory &m){
        this->_send_ewma=(u64)(m.field(0)*(1<<remycc_shift));
        this->_rec_ewma=(u64)(m.field(1)*(1<<remycc_shift));
        this->_rtt_ratio=(u64)(m.field(2)*(1<<remycc_shift));
//        this->_slow_rec_ewma=(u64)(m.field(3)*(1<<remycc_shift));
        return *this;
    }
    string str() const{
        string ret="";
        ret+=" <send,rec,rtt>=<";
        ret=ret+to_string(_send_ewma)+","\
            +to_string(_rec_ewma)+","\
            +to_string(_rtt_ratio)+"> ";
//            +to_string(_slow_rec_ewma)+"> ";
        return ret;
    };

};

struct MR{
   M _low,_upper; 
   MR():_low(),_upper(){
   };
   MR& operator=(const MemoryRange& mr){
       this-> _low=mr.getMemory(0);
       this-> _upper=mr.getMemory(1);
        return *this;
   };
   string str() const{
       return string(" <low,upper>=<")+_low.str()+","+_upper.str()+"> ";
   };
};
struct A{
    s64 _window_increment;
    u64 _window_multiple;
    u64 _intersend;
    A():_window_increment(0),_window_multiple(0),_intersend(0){};
    A& operator=(const Whisker& whisker){
        this->_window_increment=whisker.window_increment();
        this->_window_multiple=(u64)(whisker.window_multiple()*(1<<remycc_shift));
        this->_intersend=(u64)(whisker.intersend()*1000);
        return *this;
    };
    string str() const{
        string ret="";
        ret=ret+" <m,b,r>=<"+to_string(_window_multiple)+","+to_string(_window_increment)+","+to_string(_intersend)+"> ";
        return ret;
    };
};
struct WHK{
    MR _domain;
    A _action;
    int first_child;
    int child_num;
    WHK(WhiskerTree& whiskers):_domain(),_action(),first_child(0),child_num(0){
        _domain=whiskers.get_domain();
        std::vector<Whisker> whis=whiskers.leaf();
        if(whiskers.is_leaf()){
            _action=whis.front();
        }
        else{
                child_num=whiskers.num_children();
        }
    }
    string str() const{
        string ret="";
        ret+=_domain.str();
        ret+=_action.str();
        ret+="  "+to_string(first_child);
        ret+="  "+to_string(child_num);
        return ret;
    };
};



int main( int argc, char *argv[] )
{
  WhiskerTree whiskers;
  string output_filename;
  BreederOptions options;
  WhiskerImproverOptions whisker_options;
  RemyBuffers::ConfigRange input_config;
  string config_filename;

  for ( int i = 1; i < argc; i++ ) {
    string arg( argv[ i ] );
    if(arg.substr(0,3) == "hh="){
        string filename(arg.substr(3));
        int fd=open(filename.c_str(),O_RDONLY);
        if(fd<0){
            puts("err open read file");
            exit(1);
        }
        RemyBuffers::WhiskerTree tree;
        if(!tree.ParseFromFileDescriptor(fd)){
            puts("Parse err");
            exit(1);
        }
        whiskers = WhiskerTree(tree);
        filename+=".b";
        int fw=open(filename.c_str(),O_WRONLY|O_CREAT,0666);
        if(fw<0){
            puts("open err");
            exit(1);
        }
        double maxsp = 0;
        double avsp = 0;
        int cw = 0;
        int avcnt = 0;
        int id=0;
        queue<WhiskerTree> que;
        int cntnode=0;
        int cntleaf=0;
        while(!que.empty()) que.pop();
        que.push(whiskers);
        while(!que.empty()){
            WhiskerTree whtree = que.front();
            que.pop();
            cntnode++;
            WHK whk(whtree);
            int now=cntnode+que.size();
            if(whk.child_num>0)
                whk.first_child=now;
            int ssw=(int)write(fw,&whk,sizeof(whk));
            if(ssw<0)   break;
            do{
                double m = ((double)whk._action._window_multiple)/(1<<remycc_shift);;
                int b = whk._action._window_increment;
                double r = whk._action._intersend/1000.0;
                double pri = -1;
                double sp=0;
                if(whk._action._intersend){
                    sp=1.46/r;
                    avsp += sp;
                    avcnt ++;
                }
                if(whk.child_num==0 && 0<m && m<(1-eps)){
                    pri = b/(1-m);
                    printf("======id:%d,cwnd:%f\t,speed:%fMB/s\nwhisker:%s \n",id,pri, sp, whtree.str().c_str());
                }
                id++;
                if(sp>maxsp){
                    cw=pri;
                    maxsp=sp;
                }

            }while(0);
            if(whk.child_num==0) continue;
            cntleaf++;
            const std::vector<WhiskerTree> &children=whtree.children();
            for(auto x:children){
               que.push(x); 
            }
//            printf("children num:%d\n",(int)children.size());
        }
        close(fw); 
        printf("\nnode:%d,not_leaf:%d,maxsp:%f,cwnd:%d,avsp:%f\n",cntnode,cntleaf,maxsp,cw,avsp/(avcnt));
        return 0;
    }
    if ( arg.substr( 0, 3 ) == "if=" ) {
      string filename( arg.substr( 3 ) );
      int fd = open( filename.c_str(), O_RDONLY );
      if ( fd < 0 ) {
	perror( "open" );
	exit( 1 );
      }

      RemyBuffers::WhiskerTree tree;
      if ( !tree.ParseFromFileDescriptor( fd ) ) {
	fprintf( stderr, "Could not parse %s.\n", filename.c_str() );
	exit( 1 );
      }
      whiskers = WhiskerTree( tree );

      if ( close( fd ) < 0 ) {
	perror( "close" );
	exit( 1 );
      }

    } else if ( arg.substr( 0, 3 ) == "of=" ) {
      output_filename = string( arg.substr( 3 ) );

    } else if ( arg.substr( 0, 4 ) == "opt=" ) {
      whisker_options.optimize_window_increment = false;
      whisker_options.optimize_window_multiple = false;
      whisker_options.optimize_intersend = false;
      for ( char & c : arg.substr( 4 ) ) {
        if ( c == 'b' ) {
          whisker_options.optimize_window_increment = true;
        } else if ( c == 'm' ) {
          whisker_options.optimize_window_multiple = true;
        } else if ( c == 'r' ) {
          whisker_options.optimize_intersend = true;
        } else {
          fprintf( stderr, "Invalid optimize option: %c\n", c );
          exit( 1 );
        }
      }

    } else if ( arg.substr(0, 3 ) == "cf=" ) {
      config_filename = string( arg.substr( 3 ) );
      int cfd = open( config_filename.c_str(), O_RDONLY );
      if ( cfd < 0 ) {
        perror( "open config file error");
        exit( 1 );
      }
      if ( !input_config.ParseFromFileDescriptor( cfd ) ) {
        fprintf( stderr, "Could not parse input config from file %s. \n", config_filename.c_str() );
        exit ( 1 );
      }
      if ( close( cfd ) < 0 ) {
        perror( "close" );
        exit( 1 );
      }
    }
  }

  if ( config_filename.empty() ) {
    fprintf( stderr, "An input configuration protobuf must be provided via the cf= option. \n");
    fprintf( stderr, "You can generate one using './configuration'. \n");
    exit ( 1 );
  }

  options.config_range = ConfigRange( input_config );

  RatBreeder breeder( options, whisker_options );

  unsigned int run = 0;

  printf( "#######################\n" );
  printf( "Evaluator simulations will run for %d ticks\n",
    options.config_range.simulation_ticks );
  printf( "Optimizing window increment: %d, window multiple: %d, intersend: %d\n",
          whisker_options.optimize_window_increment, whisker_options.optimize_window_multiple,
          whisker_options.optimize_intersend);
  print_range( options.config_range.link_ppt, "link packets_per_ms" );
  print_range( options.config_range.rtt, "rtt_ms" );
  print_range( options.config_range.num_senders, "num_senders" );
  print_range( options.config_range.mean_on_duration, "mean_on_duration" );
  print_range( options.config_range.mean_off_duration, "mean_off_duration" );
  print_range( options.config_range.stochastic_loss_rate, "stochastic_loss_rate" );
  if ( options.config_range.buffer_size.low != numeric_limits<unsigned int>::max() ) {
    print_range( options.config_range.buffer_size, "buffer_size" );
  } else {
    printf( "Optimizing for infinitely sized buffers. \n");
  }

  printf( "Initial rules (use if=FILENAME to read from disk): %s\n", whiskers.str().c_str() );
  printf( "#######################\n" );

  if ( !output_filename.empty() ) {
    printf( "Writing to \"%s.N\".\n", output_filename.c_str() );
  } else {
    printf( "Not saving output. Use the of=FILENAME argument to save the results.\n" );
  }

  RemyBuffers::ConfigVector training_configs;
  bool written = false;

  while ( 1 ) {
    auto outcome = breeder.improve( whiskers );
    printf( "run = %u, score = %f\n", run, outcome.score );

    printf( "whiskers: %s\n", whiskers.str().c_str() );

    for ( auto &run : outcome.throughputs_delays ) {
      if ( !(written) ) {
        for ( auto &run : outcome.throughputs_delays) {
          // record the config to the protobuf
          RemyBuffers::NetConfig* net_config = training_configs.add_config();
          *net_config = run.first.DNA();
          written = true;

        }
      }
      printf( "===\nconfig: %s\n", run.first.str().c_str() );
      for ( auto &x : run.second ) {
	printf( "sender: [tp=%f, del=%f]\n", x.first / run.first.link_ppt, x.second / run.first.delay );
      }
    }

    if ( !output_filename.empty() ) {
      char of[ 128 ];
      snprintf( of, 128, "%s.%d", output_filename.c_str(), run );
      fprintf( stderr, "Writing to \"%s\"... ", of );
      int fd = open( of, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR );
      if ( fd < 0 ) {
	perror( "open" );
	exit( 1 );
      }

      auto remycc = whiskers.DNA();
      remycc.mutable_config()->CopyFrom( options.config_range.DNA() );
      remycc.mutable_optimizer()->CopyFrom( Whisker::get_optimizer().DNA() );
      remycc.mutable_configvector()->CopyFrom( training_configs );
      if ( not remycc.SerializeToFileDescriptor( fd ) ) {
	fprintf( stderr, "Could not serialize RemyCC.\n" );
	exit( 1 );
      }

      if ( close( fd ) < 0 ) {
	perror( "close" );
	exit( 1 );
      }

      fprintf( stderr, "done.\n" );
    }

    fflush( NULL );
    run++;
  }

  return 0;
}
