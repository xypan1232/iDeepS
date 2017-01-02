#include "FlagsService.h"

FlagsServiceClient::FlagsServiceClient(const std::string& id_for_flags_service) {
  this->id_for_flags_service = id_for_flags_service;
  The_FlagsService.register_flags_service_client(this);
}

FlagsServiceClient::~FlagsServiceClient() {
  for (FMap::const_iterator i=flags_traits.begin(); i!=flags_traits.end(); ++i)
    delete i->second;
  The_FlagsService.unregister_flags_service_client(this);
}
