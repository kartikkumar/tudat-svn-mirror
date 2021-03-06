/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120807    K. Kumar          File created.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_DICTIONARY_COMPARER_H
#define TUDAT_DICTIONARY_COMPARER_H

#include <algorithm>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/dictionaryEntry.h"
#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{
namespace dictionary
{

//! Dictionary comparer functor.
/*!
 * This is a functor that implements a number of overloads of the ()-operator to perform
 * various comparison operations, used in conjunction with the dictionary-based input system
 * implemented in Tudat.
 * \sa DictionaryEntry.
 */
class DictionaryComparer
{
public:

    //! Zombie constructor.
    /*!
     * Constructor that initializes all members with zombie state.
     */
    DictionaryComparer( )
        : dictionaryEntry( ),
          parameterName( "" )
    { }

    //! Constructor taking dictionary entry.
    /*!
     * Constructor that takes a dictionary entry as input and initializes the member dictionary
     * entry. The parameterName member is initialized with zombie state.
     * \param aDictionaryEntry Shared-pointer to a dictionary entry.
     * \sa DictionaryEntry.
     */
    DictionaryComparer( const input_output::dictionary::DictionaryEntryPointer aDictionaryEntry )
        : dictionaryEntry( aDictionaryEntry ),
          parameterName( "" )
    { }

    //! Constructor taking parameter name.
    /*!
     * Constructor that takes a parameter name as input and initializes the member parameter name.
     * The dictionary entry is initialized with zombie state.
     * \param aParameterName Parameter-name.
     */
    DictionaryComparer( const std::string& aParameterName )
        : dictionaryEntry( ),
          parameterName( aParameterName )
    { }

    //! Overload ()-operator to compare dictionary entries.
    /*!
     * Overloads the ()-operator to compare dictionary entries. This allows for a sorted
     * STL container of dictionary entries to be created (e.g., set).
     * \param firstDictionaryEntry Shared-pointer to a dictionary entry.
     * \param secondDictionaryEntry Shared-pointer to a dictionary entry.
     * \return True if parameter name of first dictionary entry should be sorted before parameter
     *          name of second dictionary entry.
     */
    bool operator( )( DictionaryEntryPointer firstDictionaryEntry,
                      DictionaryEntryPointer secondDictionaryEntry ) const
    {
      return firstDictionaryEntry->parameterName < secondDictionaryEntry->parameterName;
    }

    //! Overload ()-operator to compare synonym with parameter name.
    /*!
     * Overloads the ()-operator to compare a given synonym with the set parameter name. The check
     * is not case-sensitive. This allows for a list of synonyms to be checked against a specific
     * parameter name, set through the constructor.
     * \param synonym Synonym to check.
     * \return True if synonym matches parameter name (not case-sensitive).
     */
    bool operator( )( const std::string& synonym ) const
    {
        return boost::iequals( synonym, parameterName );
    }

    //! Overload ()-operator to compare dictionary entry with parameter name.
    /*!
     * Overloads the ()-operator to compare a given dictionary entry with the set parameter name.
     * The check is case-sensitive. This allows for a set of dictionary entries to be compared
     * against a specific parameter name, set through the constructor.
     * \param dictionaryEntry Dictionary entry to check.
     * \return True if dictionary entry matches parameter name (case-sensitive).
     */
    bool operator( )( DictionaryEntryPointer dictionaryEntry ) const
    {
        return !dictionaryEntry->parameterName.compare( parameterName );
    }

    //! Overload ()-operator to compare data line with dictionary entry.
    /*!
     * Overloads the ()-operator to compare a given data line with the set dictionary entry. The
     * check ensures that the paremter name and all synonyms are checked. The isCaseSensitive
     * flag in the dictionary entry is used to ensure that the comparisons are either
     * case-sensitive or not, as specified.
     * \param dataLine Data line from input stream.
     * \return True if parameter name in data line matches dictionary entry.
     */
    bool operator( )( const input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr&
                      dataLine ) const;
protected:

private:

    //! Dictionary entry.
    /*!
     * Shared-pointer to a dictionary entry.
     */
    const input_output::dictionary::DictionaryEntryPointer dictionaryEntry;

    //! Parameter name.
    /*!
     * String name of parameter.
     */
    const std::string parameterName;
};

//! Typedef for shared-pointer to DictionaryComparer object.
typedef boost::shared_ptr< DictionaryComparer > DictionaryComparerPointer;

} // namespace dictionary
} // namespace input_output
} // namespace tudat

#endif // TUDAT_DICTIONARY_COMPARER_H
