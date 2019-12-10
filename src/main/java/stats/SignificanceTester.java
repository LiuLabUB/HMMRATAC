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
import org.apache.commons.math3.distribution.*;




public class SignificanceTester {
	
	private static File input = null;
	private static File input2 = null;
	private static File output = null;
	private static int type = 0;
	private static double mu = 0.0;
	
	public static void main(String[] args) throws FileNotFoundException{
		
		for (int i = 0; i < args.length; i++) {

			switch (args[i].charAt((1))) {

			case 'i':
				input = new File(args[i + 1]);
				i++;
				break;
				
			case 'f':
				input2 = new File(args[i+1]);
				i++;
				break;
				
			case 'o':
				output = new File(args[i+1]);
				i++;
				break;
			
			case 't':
				type = Integer.parseInt(args[i+1]);
				i++;
				break;
				
			case 'm':
				mu = Integer.parseInt(args[i+1]);
				i++;
				break;
			}
		}
		double mu1 = (double)mu;
		if (input == null  || output == null || type == 0){
			printUsage();
			System.exit(0);
		}
		PrintStream OUT = new PrintStream(output);
		@SuppressWarnings("resource")
		Scanner inFile =new Scanner ((Readable) new FileReader(input));
		
		Vector<Node.MinNode> val1 = new Vector<Node.MinNode>();
		Vector<Node.MinNode> val2 = new Vector<Node.MinNode>();
		Node.MinNode temp1 = null;
		Node.MinNode temp2 = null;
		
		while (inFile.hasNext()){
			double value = inFile.nextDouble();
			temp1 = new Node.MinNode(value);
			val1.add(temp1);
		}
		
		if (type == 2){
			@SuppressWarnings("resource")
			Scanner inFile2 = new Scanner ((Readable) new FileReader(input2));
		
		
			while (inFile2.hasNext()){
				double value = inFile2.nextDouble();
				temp2 = new Node.MinNode(value);
				val2.add(temp2);
			}
		
			double[] value1 = new double[val1.size()];
			for (int i = 0;i < val1.size();i++){
				double value = val1.get(i).getMin();
				value1[i] = value;
			}
			double[] value2 = new double[val2.size()];
			for (int i = 0; i < val2.size();i++){
				double value = val2.get(i).getMin();
				value2[i] = value;
			}
		
			double pValue = gettwosidedPValue(value1,value2);
			System.out.println(pValue);
		}
		
		if (type == 1){
			double[] value1 = new double[val1.size()];
			for (int i = 0;i < val1.size();i++){
				double value = val1.get(i).getMin();
				value1[i] = value;
			}
			
			double pValue = getonesidedPValue(mu1,value1);
			System.out.println(pValue);
		}
		
		
		printTimeStamp();
	}
	
	public static void printUsage() {
		System.out.println("Usage: java -jar SignificanceTester.jar -i <InputFIleName> -f <InputFileName2> -o <OutputFileName> -t <Type> -m <mu>");
		System.out.println("Parameters:");
		System.out.println("-i <Input.txt> File containing one column of values");
		System.out.println("-f <Input2.txt> File containing one column of values");
		System.out.println("-o <Output.txt> Name of Output File");
		System.out.println("-t <Type> 1 for a one sided t-test OR 2 for a two sided t-test");
		System.out.println("-m <Mu> The value to compare for a one sided t-test");
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
