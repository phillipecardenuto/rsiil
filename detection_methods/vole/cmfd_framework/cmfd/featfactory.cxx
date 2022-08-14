/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "featfactory.h"
#include "luo.h"
#include "blur.h"
#include "dct.h"
#include "hu.h"
#include "pctfeat.h"
#include "svd.h"
#include "dwtfeat.h"
#include "fmt.h"
#include "lin.h"
#include "bravo.h"
#include "circlefeat.h"
#include "zernike.h"
#include "kpca.h"

Feature* FeatureFactory :: getInstance(featureType name, struct ft_cfg config, const BlockHandling &b, int numThreads)
{
	if (name == LUO){
		return new Luo(config, b, numThreads);
	}
	if (name == BLUR){
		return new Blur(config, b, numThreads);
	}
	if (name == DCT) {
		return new Dct(config, b, numThreads);
	}
	if (name == PCT) {
		return new Pctfeat(config, b, numThreads);
	}
	if (name == SVD){
		return new Svd(config, b, numThreads);
	}
	if (name == HU){
		return new Hu(config, b, numThreads);
	}
	if (name == DWTFEAT){
		return new Dwtfeat(config, b, numThreads);
	}
	if (name == FMT) {
		return new Fmt(config, b, numThreads);
	}
	if (name == LIN) {
		return new Lin(config, b, numThreads);
	}
	if (name == BRAVO) {
		return new Bravo(config, b, numThreads);
	}
	if (name == CIRCLE) {
		return new Circle(config, b, numThreads);
	}
	if (name == ZERNIKE) {
		return new Zernike(config, b, numThreads);
	}
	if (name == KPCA) {
		return new Kpca(config, b, numThreads);
	}

	if (name == NO || name == NOTHING){// || name == SIFT || name == SURF) {
		return new Feature(config, b, numThreads);
	}
	else {
		throw std::runtime_error("FeatureFactory::getInstance: wrong feature name specified");
	}
	return NULL;
}
