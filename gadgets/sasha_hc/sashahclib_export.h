//examplelib_export.h

#ifndef SASHAHC_EXPORT_H_
#define SASHAHC_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_SASHA_HC__)
#define EXPORTSASHAHC __declspec(dllexport)
#else
#define EXPORTSASHAHC __declspec(dllimport)
#endif
#else
#define EXPORTSASHAHC
#endif

#endif /* EXAMPLE_EXPORT_H_ */