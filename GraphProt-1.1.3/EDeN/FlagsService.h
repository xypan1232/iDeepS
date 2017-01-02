/* -*- mode:c++ -*-
 * FlagsService.h
 *
 *  Created on: Aug 19, 2010
 *      Author: Paolo Frasconi
 */

#ifndef FLAGSSERVICE_H_
#define FLAGSSERVICE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <stdexcept>
#include <map>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <typeinfo>
#include <vector>

/** 
 * @class FlagsTraits
 * 
 * @brief Contains traits of flags. The name of a flag is used by
 * Prolog's predicates set_klog_flag and get_klog_flag. The
 * description is for the interactive help.
 */
class FlagsTraits {
private:
  std::string name; 
  std::string description;
public:
  FlagsTraits(const std::string& name, const std::string& description) {
    this->name = name;
    this->description = description;
  }
  std::string& get_name(void) { return name; }
  std::string& get_description(void) { return description; }
  virtual void dummy(void) {}; // make dynamic_cast possible...
};
class BoolFlagsTraits : public FlagsTraits {
private:
  bool* ref;
public:
  BoolFlagsTraits(bool* ref, const std::string& name, const std::string& description):
    FlagsTraits(name,description) {
    this->ref = ref;
  }
  bool* get_ref(void) {return ref;}
};
class IntFlagsTraits : public FlagsTraits {
private:
  int* ref;
public:
  IntFlagsTraits(int* ref, const std::string& name, const std::string& description):
    FlagsTraits(name,description) {
    this->ref = ref;
  }
  int* get_ref(void) {return ref;}
};
class UnsignedFlagsTraits : public FlagsTraits {
private:
  unsigned* ref;
public:
  UnsignedFlagsTraits(unsigned* ref, const std::string& name, const std::string& description):
    FlagsTraits(name,description) {
    this->ref = ref;
  }
  unsigned* get_ref(void) {return ref;}
};
class DoubleFlagsTraits : public FlagsTraits {
private:
  double* ref;
public:
  DoubleFlagsTraits(double* ref, const std::string& name, const std::string& description):
    FlagsTraits(name,description) {
    this->ref = ref;
  }
  double* get_ref(void) {return ref;}
};
class StringFlagsTraits : public FlagsTraits {
private:
  std::string* ref;
public:
  StringFlagsTraits(std::string* ref, const std::string& name, const std::string& description):
    FlagsTraits(name,description) {
    this->ref = ref;
  }
  std::string* get_ref(void) {return ref;}
};

/**
 * @class FlagsServiceClient
 *
 * Interface for clients of the Flags service. Every class that needs
 * to use the Flags service should inherit from this. A flag is a
 * prolog reference to a data member of the class using the
 * service. The class is expected to push its own flags using the
 * method new_flag. Currently only bool, int, unsigned, double, and
 * std::string are valid types for data members accessible via
 * flags. In their constructors, derived classes should initialize the
 * base class with a string (id) that is used to identify the specific
 * object of that class. Then Prolog can communicate with the class
 * via the predicates set_klog_flag(+Id,+FlagName,+FlagValue) and
 * get_klog_flag(+Id,+FlagName,-FlagValue). Prolog values are cast
 * into the right type of the object variable automatically. The
 * derived class is responsible for raising exceptions if the values
 * are not as expected. Exceptions for invalid flag names are raised
 * here.
 */
class FlagsServiceClient {
private:
  std::string id_for_flags_service;
  typedef std::map<std::string,FlagsTraits*> FMap;
  FMap flags_traits;
public:
  std::string id_str(void) const {return id_for_flags_service;}
  void new_flag(bool* ref, const std::string& name, const std::string& description) {
    flags_traits[name] = new BoolFlagsTraits(ref,name,description);
  }
  void new_flag(int* ref, const std::string& name, const std::string& description) {
    flags_traits[name] = new IntFlagsTraits(ref,name,description);
  }
  void new_flag(unsigned* ref, const std::string& name, const std::string& description) {
    flags_traits[name] = new UnsignedFlagsTraits(ref,name,description);
  }
  void new_flag(double* ref, const std::string& name, const std::string& description) {
    flags_traits[name] = new DoubleFlagsTraits(ref,name,description);
  }
  void new_flag(std::string* ref, const std::string& name, const std::string& description) {
    flags_traits[name] = new StringFlagsTraits(ref,name,description);
  }
  std::string document_flags(void) const {
    std::ostringstream outs;
    outs.setf (std::ios::left, std::ios::adjustfield);
    for (FMap::const_iterator i=flags_traits.begin(); i!=flags_traits.end(); ++i) {
      outs << std::setw(23) << id_for_flags_service + "::" + i->second->get_name() << " [";
      BoolFlagsTraits* pb; IntFlagsTraits* pi; UnsignedFlagsTraits* pu;  DoubleFlagsTraits* pd; StringFlagsTraits* ps;
      if ((pb=dynamic_cast<BoolFlagsTraits*>(i->second)) != NULL)
        outs << *(pb->get_ref());
      else if ((pi=dynamic_cast<IntFlagsTraits*>(i->second)) != NULL)
        outs << *(pi->get_ref());
      else if ((pu=dynamic_cast<UnsignedFlagsTraits*>(i->second)) != NULL)
        outs << *(pu->get_ref());
      else if ((pd=dynamic_cast<DoubleFlagsTraits*>(i->second)) != NULL)
        outs << *(pd->get_ref());
      else if ((ps=dynamic_cast<StringFlagsTraits*>(i->second)) != NULL)
        outs << *(ps->get_ref());
      else
        throw std::runtime_error("Nasty bug - unrecognized type in document_flags");
      outs << "]" << std::endl << i->second->get_description() << std::endl << std::endl;
    }
    return outs.str();
  }
  virtual ~FlagsServiceClient();
  // virtual void set_flag_callback(const std::string& flag_name, const std::string& flag_value) = 0;
  // virtual std::string get_flag_callback(const std::string& flag_name) const = 0;

  void set_flag(const std::string& flag_name, const std::string& flag_value) {
    FMap::iterator found = flags_traits.find(flag_name);
    if (found==flags_traits.end())
      throw std::runtime_error("Unknown flag name");
    BoolFlagsTraits* pb; IntFlagsTraits* pi; UnsignedFlagsTraits* pu;  DoubleFlagsTraits* pd; StringFlagsTraits* ps;
    if ((pb=dynamic_cast<BoolFlagsTraits*>(found->second)) != NULL)
      *(pb->get_ref()) = (flag_value=="true"||flag_value=="on"||flag_value=="yes");
    else if ((pi=dynamic_cast<IntFlagsTraits*>(found->second)) != NULL)
      *(pi->get_ref()) = atoi(flag_value.c_str());
    else if ((pu=dynamic_cast<UnsignedFlagsTraits*>(found->second)) != NULL)
      *(pu->get_ref()) = unsigned(atoi(flag_value.c_str()));
    else if ((pd=dynamic_cast<DoubleFlagsTraits*>(found->second)) != NULL)
      *(pd->get_ref()) = atof(flag_value.c_str());
    else if ((ps=dynamic_cast<StringFlagsTraits*>(found->second)) != NULL)
      *(ps->get_ref()) = flag_value;
    else
      throw std::runtime_error("Nasty bug - unrecognized type in set_flag");
  }
  std::string get_flag(const std::string& flag_name) {
    FMap::iterator found = flags_traits.find(flag_name);
    std::stringstream oss;
    if (found==flags_traits.end())
      throw std::runtime_error("Unknown flag name");
    BoolFlagsTraits* pb; IntFlagsTraits* pi; UnsignedFlagsTraits* pu;  DoubleFlagsTraits* pd; StringFlagsTraits* ps;
    if ((pb=dynamic_cast<BoolFlagsTraits*>(found->second)) != NULL)
      oss << *(pb->get_ref());
    else if ((pi=dynamic_cast<IntFlagsTraits*>(found->second)) != NULL)
      oss << *(pi->get_ref());
    else if ((pu=dynamic_cast<UnsignedFlagsTraits*>(found->second)) != NULL)
      oss << *(pu->get_ref());
    else if ((pd=dynamic_cast<DoubleFlagsTraits*>(found->second)) != NULL)
      oss << *(pd->get_ref());
    else if ((ps=dynamic_cast<StringFlagsTraits*>(found->second)) != NULL)
      oss << *(ps->get_ref());
    else
      throw std::runtime_error("Nasty bug - unrecognized type in get_flag");
    return oss.str();
  }
  FlagsServiceClient(const std::string& id_for_flags_service);
};

/**
 * @class FlagsService
 *
 * A singleton storing all the clients of the flags service. Prolog
 * calls for seting, gettin, and documenting flags are intercepted
 * here and dispatched to the right client by looking up the client's
 * id in the internal storage of this object.
 */

class FlagsService {

public:
  typedef std::map<std::string,FlagsServiceClient*> Map;
private:

  Map pool;

  FlagsService() {};
  FlagsService(FlagsService const& rhs);
  FlagsService& operator=(FlagsService const& rhs);
  FlagsServiceClient* find_range_checked(const std::string& client_id) const {
    Map::const_iterator found = pool.find(client_id);
    if ( found == pool.end() )
      throw std::out_of_range("Flags client identifier not registered");
    return found->second;
  }

public:
  static FlagsService& get_instance() {
	  static FlagsService instance;
	  return instance;
  }

  // static unsigned to_unsigned(const std::string& s) { return unsigned(atoi(s.c_str())); }
  // static int to_int(const std::string& s) { return atoi(s.c_str()); }
  // static double to_double(const std::string& s) { return atof(s.c_str()); }
  // static bool to_bool(const std::string& s) { 
  //   std::string r=s;
  //   std::transform(r.begin(), r.end(), r.begin(), ::tolower);
  //   return (r=="true" || r=="yes" || r=="on" || r=="t");
  // }
  // template <typename T>
  // static std::string to_string(const T x) {
  //   std::stringstream oss; oss << x; return oss.str();
  // }

  // Register a new client in the pool. The client must respond to the
  // FlagsClient interface defined above.
  void register_flags_service_client(FlagsServiceClient* c) {
    Map::iterator found = pool.find(c->id_str());
    if ( found != pool.end() ) {
      throw std::out_of_range("Duplicate client identifier in register_client()");
    }
    pool[c->id_str()] = c;
  }
  void unregister_flags_service_client(FlagsServiceClient* c) {
    Map::iterator found = pool.find(c->id_str());
    if ( found == pool.end() ) {
      throw std::out_of_range("Unknown client identifier in unregister_client()");
    }
    pool.erase(c->id_str());
  }

  void set_flag(const std::string client_id, const std::string& flag_name, const std::string& flag_value) {
    //find_range_checked(client_id)->set_flag_callback(flag_name, flag_value);
    find_range_checked(client_id)->set_flag(flag_name, flag_value);
  }
  std::string get_flag(const std::string client_id, const std::string& flag_name) {
    // return find_range_checked(client_id)->get_flag_callback(flag_name);
    return find_range_checked(client_id)->get_flag(flag_name);
  }
  std::string document_flags(const std::string client_id) {
    return find_range_checked(client_id)->document_flags();
  }

  bool is_registered_client(const std::string& client_id) const {
    Map::const_iterator found = pool.find(client_id);
    return ( found != pool.end() );
  }

  void get_client_ids(std::vector<std::string>& ids) const {
    for(Map::const_iterator it = pool.begin(); it != pool.end(); ++it) {
      ids.push_back(it->first);
    }
  }

  void get_client_ids_alt(Map::const_iterator& begin, Map::const_iterator& end) const {
	begin = pool.begin();
	end = pool.end();
  }

};

extern FlagsService& The_FlagsService;



#endif /* FLAGSSERVICE_H_ */
