#ifndef __RUNDBINTERFACE__HH__
#define __RUNDBINTERFACE__HH__

// BE SURE TO INCLUDE <mysql.h> IN YOUR SEARCH PATH
// FOR ACLIC, EXECUTE ".include /usr/include/mysql"

#include "TMySQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TMap.h"
#include "TObjString.h"
#include <iostream>

class RunDBInterface {

 public:

  RunDBInterface(TString address, TString database, TString table, TString user, TString password);
  void LoadRun(Int_t run_number);
  TString Get(TString key);

//  // EXAMPLE:
//  RunDBInterface *my = new RunDBInterface("mysqlserver","rundb_v2_analysis","cef3","cmsdaq_ro","....");
//  my->LoadRun(462);
//  my->GetMap()->Print();
//  cout << "WARNING: ALL METHODS RETURN TSTRINGS, CONVERT THEM IF APPROPRIATE" << endl;
//  cout << my->Get("run_nevents") << end;

  ~RunDBInterface();
  TSQLServer* GetDBObject();
  TMap* GetMap();

 private:

  TSQLServer *db;
  TSQLResult *res;
  TMap *mymap;
  TString table_;

};

#endif

