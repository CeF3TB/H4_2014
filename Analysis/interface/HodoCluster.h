#ifndef HodoCluster_h
#define HodoCluster_h

class HodoCluster {

 public:

  HodoCluster( int totFibres, float width ) {
    size_ = 0;
    pos_ = 0.;
    hodoTotFibres_ = totFibres;
    fibreWidth_ = width;
  }

  ~HodoCluster() {};

  int getSize() { return size_; };
  float getPosition() { return pos_; };

  void addFibre( int i );


 private:

  int size_;
  float pos_;

  int hodoTotFibres_;
  float fibreWidth_; // in mm

};


#endif
