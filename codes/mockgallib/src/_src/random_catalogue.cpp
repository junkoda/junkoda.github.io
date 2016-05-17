#include "random_catalogue.h"

static void fill_random_satellites(Snapshot const * const snp,
			   Sky const * const sky,
			   Remap const * const remap,
			   Slice const * const slice,
			   Catallogue* const cat);

void random_catalogue_generate_satellites(Hod* const hod,
					  Snapshots const * const snapshots,
					  Sky const * const sky,
					  Catalogue * const cat)
{
  for(Snapshots::const_iterator snp= snapshots->begin();
      snp != snapshots->end(); ++snp) {
    
    fill_random_satellites(*snp, sky, slice, cat);
  }
}

void fill_random_satellites(Hod* const hod,
			    Sky const * const sky,
			    Catalogue* const cat)
{
  
}
