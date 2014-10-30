#ifndef __RUNDBINTERFACE__CC__
#define __RUNDBINTERFACE__CC__

#include "RunDBInterface.h"
#include <iostream>

using namespace std;

RunDBInterface::RunDBInterface(TString address, TString database, TString table, TString user, TString password){
  mymap = new TMap();
  db = TSQLServer::Connect(Form("mysql://%s/%s",address.Data(),database.Data()), user.Data(), password.Data());
  if (!db) {
    cout << "Error in opening MySQL database" << endl;
    return;
  }
  table_ = table;
}

RunDBInterface::~RunDBInterface(){
  db->Close();
  if (mymap) delete mymap;
  if (res) delete res;
  if (db) delete db;
}

TSQLServer* RunDBInterface::GetDBObject(){
  return db;
}

void RunDBInterface::LoadRun(Int_t run_number){
  mymap->Clear();
  res = db->Query(Form("SELECT * FROM %s WHERE run_number=%d",table_.Data(),run_number));
  if (!res) {
    cout << "Impossible to perform query on table " << table_.Data() << endl;
    return;
  }
  int nfields = res->GetFieldCount();
  int nrows = res->GetRowCount();

  for (int i=0; i < nrows; i++) {
    TSQLRow *row = res->Next();
    for (int j = 0; j < nfields; j++) {
      TObjString *key = new TObjString(res->GetFieldName(j));
      TObjString *val = new TObjString(row->GetField(j));
      mymap->Add((TObject*)key,(TObject*)val);
    }
  }
}

TMap* RunDBInterface::GetMap(){
  return mymap;
}

TString RunDBInterface::Get(TString key){
  TObjString* val = (TObjString*)(GetMap()->GetValue(key.Data()));
  if (!val) {
    cout << "key not found: " << key.Data() << endl;
    return TString();
  }
  return val->String();
}


#endif
