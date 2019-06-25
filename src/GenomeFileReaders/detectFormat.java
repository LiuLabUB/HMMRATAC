package GenomeFileReaders;

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

public class detectFormat {
	
	private String fileName;
	
	
	public detectFormat(String f){
		fileName = f.toLowerCase();
	}
	
	public boolean detectSAM(){
		if (fileName.endsWith("sam")){return true;}
		else{return false;}
	}
	public boolean detectBAM(){
		if (fileName.endsWith("bam")){return true;}
		else{return false;}
	}
	public boolean detectBED(){
		if (fileName.endsWith("bed")){return true;}
		else {return false;}
	}
	
	
}
