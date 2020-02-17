 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <string>
#include "ThirdPartyHeadersEnd.h"
#include "SzlFileLoader.h"
#include "stringformat.h"
namespace tecplot { namespace ___3933 { class IODescription { public: typedef uint32_t SegmentIndex_t; static ___4352 const               NO_VAR = BAD_VAR_INDEX; static ___4636 const              NO_ZONE = BAD_ZONE_INDEX; static ___2090::___2980 const NO_PARTITION = ___2090::INVALID_PARTITION; static SegmentIndex_t const           NO_SEGMENT = SegmentIndex_t(-1); private: char const*              ___2495; ___4352               m_var; ___4636              ___2677; ___2090::___2980 m_partition; SegmentIndex_t           m_segment; char const*              m_suffix; public: explicit IODescription( char const*              ___2685 = NULL, ___4352               ___4336 = NO_VAR, ___4636              zone = NO_ZONE, ___2090::___2980 ___2977 = NO_PARTITION, SegmentIndex_t           segment = NO_SEGMENT, char const*              suffix = NULL) : ___2495(___2685) , m_var(___4336) , ___2677(zone) , m_partition(___2977) , m_segment(segment) , m_suffix(suffix) {} char const* ___2685() const { return ___2495; } ___4352 ___4336() const { return m_var; } ___4636 zone() const { return ___2677; } ___2090::___2980 ___2977() const { return m_partition; } SegmentIndex_t segment() const { return m_segment; } char const* suffix() const { return m_suffix; } ___372 ___2067() const { return ___4226; } ___372 isEmpty() const { return ___2495==NULL; } void getFormattedDescription( char*  formattedDescription, size_t formattedDescriptionSize) const { ___372 isAsciiOnly = ___1305; size_t ___2865; if ( ___2495 != NULL ) { ___2865 = snprintf(formattedDescription, formattedDescriptionSize, "%s", ___2495); if ( ___2865 > 0 && formattedDescription[___2865-1] == '*' ) { ___2865--; formattedDescription[___2865] = '\0'; isAsciiOnly = ___4226; } } else ___2865 = snprintf(formattedDescription, formattedDescriptionSize, "unspecified"); if ( m_var != NO_VAR && ___2865 < formattedDescriptionSize ) ___2865 += snprintf(formattedDescription+___2865, formattedDescriptionSize-___2865, "%sVar%" PRIu64, ___2495 != NULL ? "For" : "", uint64_t(m_var+1)); if ( ___2677 != NO_ZONE && ___2865 < formattedDescriptionSize ) ___2865 += snprintf(formattedDescription+___2865, formattedDescriptionSize-___2865, "%sZone%" PRIu64, ___2495 != NULL && m_var == NO_VAR ? "For" : "", uint64_t(___2677+1)); if ( m_partition != NO_PARTITION && ___2865 < formattedDescriptionSize ) ___2865 += snprintf(formattedDescription+___2865, formattedDescriptionSize-___2865, "%sPartition%" PRIu64, ___2495 != NULL && m_var == NO_VAR && ___2677 == NO_ZONE ? "For" : "", uint64_t(m_partition+1)); if ( m_segment != NO_SEGMENT && ___2865 < formattedDescriptionSize ) ___2865 += snprintf(formattedDescription+___2865, formattedDescriptionSize-___2865, "%sSegment%" PRIu64, ___2495 != NULL && m_var == NO_VAR && ___2677 == NO_ZONE && m_partition == NO_PARTITION ? "For" : "", uint64_t(m_segment+1)); if ( m_suffix != NULL && ___2865 < formattedDescriptionSize ) ___2865 += snprintf(formattedDescription+___2865, formattedDescriptionSize-___2865, "%s", m_suffix); if ( isAsciiOnly && ___2865 < formattedDescriptionSize ) { formattedDescription[___2865] = '*'; if ( ___2865 < formattedDescriptionSize ) ___2865++; formattedDescription[___2865] = '\0'; } } }; }}
