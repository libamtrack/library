package amtrackservice.server;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * data file used as input/output for libamtrack wrapper
 * @author Christoph
 * 
 */
public class CalculationFile {

	private static PrintWriter writer;

	/**
	 * read data from file and translate them to hashmap
	 * @param path
	 * @return 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static HashMap<String, String> readData(String path)
			throws FileNotFoundException, IOException {
		HashMap<String, String> keyVal = readStrings(path);
		return keyVal;
	}
	
	/**
	 * read data from file and return them as string
	 * @param fileName
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static String readContent(String fileName)
	throws FileNotFoundException, IOException {

		String result = "";
		BufferedReader reader = null;

		try {
			reader = new BufferedReader(new FileReader(fileName));
		} catch (IOException e) {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException ex) {
			}
			throw new FileNotFoundException("File " + fileName + " not found");
		}

		try {
			String nextLine = reader.readLine();
			while (nextLine != null) {
				result += nextLine + "\n";
				nextLine = reader.readLine();
			}

		} catch (FileNotFoundException e) {
			throw e;
		} catch (IOException e) {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException ex) {
				}
			}
			throw e;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException ex) {
				}
			}
		}

		return result;
	}

	/**
	 * Write data to the file
	 * 
	 * @param datei 
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public static void writeData(String file, String path) throws IOException,
			FileNotFoundException {

		createFile(path);
		String[] split = file.split("\n");
		for(String s:split){
			writeLine(s);
		}
		closeFile();
	}

	/**
	 * delete file
	 * @param path
	 */
	public static void deleteFile(String path) {
		File file = new File(path);
		if (!file.delete()) {
			if (file.isDirectory()) {
				File files[] = file.listFiles();
				for (int i = 0; i < files.length; i++)
					deleteFile(files[i].getAbsolutePath());
				file.delete();
			}
		}
	}

	/**
	 * TODO
	 * @param path
	 * @return ArrayList<String> 
	 */
	public static ArrayList<String> readOrder(String path) {
		File file = new File(path);
		ArrayList<String> dataVector = new ArrayList<String>();
		if (file.isDirectory()) {
			File[] files = file.listFiles();
			for (int i = 0; i < files.length; i++) {
				try {
					if(files[i].isDirectory())
					dataVector.add(files[i].getCanonicalPath());
				} catch (IOException e) {
					System.out.println("Data cannot be read");
				}
			}
		}
		return dataVector;
	}

	/**
	 * Closes the currently open file.
	 */
	private static void closeFile() throws IOException {
		if (writer == null)
			throw new IOException("Cannot close file, file not open.");
		writer.flush();
		writer.close();
		writer = null;
	}

	/**
	 * Either creates a new file or overwrites the existing File with a empty
	 * file, and opens the new File for Data input. Data can now be written to
	 * the file via writeLine(). The File must then be closed with closeFile();
	 */
	private static void createFile(String fileName) throws IOException,
			FileNotFoundException {
		File file = new File(fileName);
		File dir = file;
		if (dir.isFile())
			dir = file.getParentFile();
		dir.mkdirs();
		writer = null;

		// if File does already exist, overwrite
		if (file.exists())
			file.delete();

		// try to create File
		try {
			file.createNewFile();
		} catch (IOException e1) {
			throw new IOException("Cannot create file " + fileName);
		}

		// Cancel if File cannot be written to.
		if (!file.canWrite()) {
			throw new IOException("Cannot write to file " + fileName);
		}

		// create Stream:
		try {
			writer = new PrintWriter(new BufferedWriter(new FileWriter(file)));
		} catch (FileNotFoundException e) {
			throw e;

		} catch (IOException e) {
			throw e;
		}

	}

	/**
	 * Read data from file to HashMap
	 * 
	 * @return 
	 */
	private static HashMap<String, String> readStrings(String fileName)
			throws FileNotFoundException, IOException {

		HashMap<String, String> result = new HashMap<String, String>();
		BufferedReader reader = null;

		try {
			reader = new BufferedReader(new FileReader(fileName));
		} catch (IOException e) {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException ex) {
			}
			throw new FileNotFoundException("File " + fileName + " not found");
		}

		try {
			String nextLine = reader.readLine();
			while (nextLine != null) {

				String[] keyVal;
				keyVal = nextLine.split(":");
				if (keyVal.length == 2) {
					result.put(keyVal[0].trim(), keyVal[1].trim());
				}
				nextLine = reader.readLine();
			}

		} catch (FileNotFoundException e) {
			throw e;
		} catch (IOException e) {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException ex) {
				}
			}
			throw e;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException ex) {
				}
			}
		}

		return result;
	}

	/**
	 * Writes a line to the file, if a file is currently opened for writing
	 * 
	 * @param line 	The Line to be written to the File
	 */
	private static void writeLine(String line) throws IOException {
		if (writer == null)
			throw new IOException("Cannot write to file, file not open");
		writer.println(line);
	}

}