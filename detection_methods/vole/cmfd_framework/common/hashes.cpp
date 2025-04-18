/*
	Copyright(c) 2012 Christian Riess <christian.riess@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#include "hashes.h"

namespace vole {

unsigned long Hashes::getHash(const char *str, HashMethod m) {
	if (m == HASH_djb2) return Hashes::djb2(str);
	if (m == HASH_sdbm) return Hashes::sdbm(str);
	return 0;
}

unsigned long Hashes::djb2(const char *str)
{
	unsigned long hash = 5381;
	int c;
	c = *str++;
	while (c != 0) {
		if (c < 0) c += 256;
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
		c = *str++;
	}

	return hash;
}


unsigned long Hashes::sdbm(const char *str)
{
	unsigned long hash = 0;
	int c;
	c = *str++;
	while (c != 0) {
		if (c < 0) c += 256;
		hash = c + (hash << 6) + (hash << 16) - hash;
		c = *str++;
	}

	return hash;
}

}

