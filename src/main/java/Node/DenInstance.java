package Node;
/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

import net.sf.javaml.core.DenseInstance;

public class DenInstance extends DenseInstance {
	
	private String _id;
	
	public DenInstance(double[] att, String ID) {
		super(att);
		_id = ID;
	}
	
	
	public String getStringID() {
		
		return _id;
	}
	
}
