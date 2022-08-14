/*
	Copyright(c) 2012 Johannes Jordan <johannes.jordan@cs.fau.de>.

	This file may be licensed under the terms of of the GNU General Public
	License, version 3, as published by the Free Software Foundation. You can
	find it here: http://www.gnu.org/licenses/gpl.html
*/

#ifndef PROGRESS_OBSERVER_H
#define PROGRESS_OBSERVER_H

namespace vole {

	class ProgressObserver
	{
	public:
		/// if false: cancel job
		virtual bool update(int percent) = 0;
	};
}

#endif // PROGRESS_OBSERVER_H
