package stats;
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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.Scanner;
import java.util.Vector;

import org.apache.commons.math3.stat.inference.TestUtils;



public class ColumnCounter{
	public static File input1 = null;
	public static File input2 = null;
	public static File OutPut = null;
	public static String type = "";
	
	
	public static void main(String[] args) throws FileNotFoundException{
	
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {

			case 'i':
				input1 = new File(args[i + 1]);
				i++;
				break;
				
			case 'f':
				input2 = new File(args[i+1]);
				i++;
				break;
			
			case 'o':
				OutPut = new File(args[i+1]);
				i++;
				break;
				
			case 't':
				type = new String(args[i+1]);
				i++;
				break;
				
			
			
			}
		}
		
		if (input1 == null || input2 == null || OutPut == null || type == ""){
			printUsage();
		}
		
		
		//needs to be able to look at multiple columns and compare them. not just one to two but also two to three etc
		//also - columns are sometimes of different lengths.  There WILL HAVE to be TWO inputfiles
		PrintStream OUT = new PrintStream(OutPut);
		@SuppressWarnings("resource")
		Scanner inFile =new Scanner ((Readable) new FileReader(input1));
		@SuppressWarnings("resource")
		Scanner inFile2 = new Scanner ((Readable) new FileReader(input2));
		
		String word = "word";
		String num = "number";
		
		if (type.equals(word)){
			Vector<String> Val1 = new Vector<String>();
			Vector<String> Val2 = new Vector<String>();
		
			while (inFile.hasNext()){
				
				String value1 = inFile.next();
				//String value2 = inFile.next();
				//System.out.println(value1+"\t"+value2);
				Val1.add(value1);
				//Val2.add(value2);
				
			}
			while (inFile2.hasNext()){
				String value2 = inFile2.next();
				Val2.add(value2);
				
			}
			
			for (int i = 0;i < Val1.size();i++){
				String value1 = Val1.get(i);
				int counter = 0;
				for (int a = 0;a < Val2.size();a++){
					String value2 = Val2.get(a);
					if (value1.equals(value2)){
						counter++;
					}
				}
				OUT.println(value1+"\t"+counter);
			
			}
			
		}
		if (type.equals(num)){
			Vector<Node.MinNode> Val1 = new Vector<Node.MinNode>();
			Node.MinNode temp1 = null;
			Vector<Node.MinNode> Val2 = new Vector<Node.MinNode>();
			Node.MinNode temp2 = null;
			
			while (inFile.hasNext()){
				double value1 = inFile.nextDouble();
				//double value2 = inFile.nextDouble();
				temp1 = new Node.MinNode(value1);
				//temp2 = new Node.MinNode(value2);
				Val1.add(temp1);
				//Val2.add(temp2);
			}
			while (inFile2.hasNext()){
				double value2 = inFile2.nextDouble();
				temp2 = new Node.MinNode(value2);
				Val2.add(temp2);
			}
			for (int i = 0;i < Val1.size();i++){
				double value1 = Val1.get(i).getMin();
				int counter = 0;
				for (int a = 0;a < Val2.size();a++){
					double value2 = Val2.get(a).getMin();
					if (value1 == value2){
						counter++;
					}
				}
				OUT.println(value1+"\t"+counter);
			
			}
			
			
		}
		
		
		
		
		
		printTimeStamp();
	}
	
	public static void printUsage() {
		System.out.println("Usage: java -jar ColumnCounter.jar -i <InputFIleName> -f <InputFileName2> -o <OutputFileName> -t <Type>");
		System.out.println("Parameters:");
		System.out.println("-i <Input.txt> File containing one column of values");
		System.out.println("-f <Input2.txt> File containing one column of values");
		System.out.println("-o <Output.txt> Name of Output File");
		System.out.println("-t <Type> Can be either <word> or <number>");
	}
	
	private static double getonesidedPValue(double mu,double[] sample1){
		double PValue = TestUtils.tTest(mu, sample1);
		
		return PValue;
	}
	
	private static double gettwosidedPValue(double[] sample1,double[] sample2){
		
		double PValue = TestUtils.tTest(sample1, sample2);
		return PValue;
		
	
	}
	
	private static void printTimeStamp() {
		//Outputs a time stamp onto Standard Output
		Calendar timestamp = new GregorianCalendar();
		System.out.println("Program Complete");
		System.out.print("Current Time: " + (timestamp.get(Calendar.MONTH)+1) + "-" + timestamp.get(Calendar.DAY_OF_MONTH) + "-" + timestamp.get(Calendar.YEAR));
		System.out.println("\t" + timestamp.get(Calendar.HOUR_OF_DAY) + ":" + timestamp.get(Calendar.MINUTE) + ":" + timestamp.get(Calendar.SECOND));
	}
}
