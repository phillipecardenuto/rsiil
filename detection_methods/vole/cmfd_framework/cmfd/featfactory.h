/*
	Copyright(c) 2011 Vincent Christlein <vincent.christlein@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef FEATFACTORY_H
#define FEATFACTORY_H

#include <string>

#include "copymoveframework.h"
#include "feature.h"

class FeatureFactory {
	public:
        FeatureFactory(){}
        ~FeatureFactory(){}
		/// returns the proper Featureinstance
		static Feature* getInstance(featureType name, struct ft_cfg, const BlockHandling &b, int numThreads = -1);
};
#endif // FEATFACTORY_H
