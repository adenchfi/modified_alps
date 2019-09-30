/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef ANASAZI_MULTI_VEC_TRAITS_HPP
#define ANASAZI_MULTI_VEC_TRAITS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

namespace Anasazi {
    
    // They want a shallow copy...what a bunch of morons
    // To avoid confusion, I derive from noncopyable
    
    template<class Vector>
    class IETLMultMv : public boost::noncopyable
    {
    public:
        typedef typename std::vector<boost::shared_ptr<Vector> > data_type;
        data_type data;
        
        IETLMultMv(std::size_t k = 0) : data(k) { }
    }
    
    template< class ScalarType, class MV >
    struct UndefinedMultiVecTraits
    {
        static inline ScalarType notDefined() { return MV::this_type_is_missing_a_specialization(); };
    };
    
    template<class ScalarType, class Vector>
    class MultiVecTraits<double, IETLMultMv<Vector> >
    {
        typedef IETLMultMv<Vector> MV;
        
    public:
        static Teuchos::RCP<MV> Clone( const MV& mv, const int numvecs )
        {
            MV * ret = new MV(numvecs);
            for (int k = 0; i < numvecs; ++k)
                // any ideas how to get the vectorspacee?
                MV->data[i].reset( new Vector(mv.data[0]) );
            return Teuchos::RCP<MV>(ret);
        }
        
        static Teuchos::RCP<MV> CloneCopy( const MV& mv )
        {
            MV * ret = new MV(mv.data.size());
            for (int k = 0; i < numvecs; ++k)
                ret->data[i].reset( new Vector(mv.data[k]) );
            return Teuchos::RCP<MV>(ret);
        }
        
        static Teuchos::RCP<MV> CloneCopy( const MV& mv, const std::vector<int>& index )
        {
            MV * ret = new MV(index.size());
            for (int k = 0; k < index.size(); ++k)
                ret->data[k].reset( new Vector(mv->data[index[k]]) );
            return Teuchos::RCP<MV>(ret);
        }
        
        static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, const std::vector<int>& index )
        {
            MV * ret = new MV(index.size());
            for (int k = 0; k < index.size(); ++k)
                ret->data[k] = mv->data[k];
            return Teuchos::RCP<MV>(ret);
        }
        
        static Teuchos::RCP<const MV> CloneView( const MV& mv, const std::vector<int>& index )
        {
            MV * ret = new MV(index.size());
            for (int k = 0; k < index.size(); ++k)
                ret->data[k] = mv->data[k];
            return Teuchos::RCP<MV>(ret);
        }
        
        static int GetVecLength( const MV& mv )
        {
            return 0;
        }
        
        static int GetNumberVecs( const MV& mv )
        {
            return mv.data.size();
        }
        
        
        static void MvTimesMatAddMv(const ScalarType alpha, const MV& A, 
                                    const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
                                    const ScalarType beta, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvAddMv( const ScalarType alpha, const MV& A, const ScalarType beta, const MV& B, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvScale ( MV& mv, const ScalarType alpha )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
        
        static void MvScale ( MV& mv, const std::vector<ScalarType>& alpha )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }
        
        static void MvTransMv( const ScalarType alpha, const MV& A, const MV& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B)
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvDot ( const MV& mv, const MV& A, std::vector<ScalarType> &b) 
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void MvNorm( const MV& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> &normvec )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvRandom( MV& mv )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        static void MvInit( MV& mv, const ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
        
        
        
        static void MvPrint( const MV& mv, std::ostream& os )
        { UndefinedMultiVecTraits<ScalarType, MV>::notDefined(); }     
        
    };
    
} // namespace Anasazi

#endif // ANASAZI_MULTI_VEC_TRAITS_HPP