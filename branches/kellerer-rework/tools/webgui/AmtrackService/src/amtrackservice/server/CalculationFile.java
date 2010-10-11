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
 * Liest aus einer bestehenden txt-Datei und schreibt in einer txt-Datei liest
 * Ordner aus
 * 
 * @author Christoph
 * 
 */
public class CalculationFile {

	private static PrintWriter writer;

	/**
	 * Liest die Datei und beschreibt ein Untersuchung Objekt mit den Daten
	 * 
	 * @param pfad
	 * 
	 * @return das aus der Datei ausgelesene Untersuchung Objekt
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static HashMap<String, String> dateiAuslesen(String pfad)
			throws FileNotFoundException, IOException {
		HashMap<String, String> keyVal = readStrings(pfad);

		return keyVal;
	}

	/**
	 * Schreibt die Daten des übergebenen Speicherbaren Objekts in eine Datei
	 * 
	 * @param datei
	 *            Das Speicherbare Objekt.
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public static void dateiSchreiben(String file, String path) throws IOException,
			FileNotFoundException {

		createFile(path);
		String[] split = file.split("\n");
		for(String s:split){
			writeLine(s);
		}
		closeFile();
	}

	/**
	 * Falls die Datei existiert, wird diese gelöscht.
	 * 
	 * @param pfad
	 *            der Datei
	 */
	public static void deleteFile(String pfad) {
		File file = new File(pfad);

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
	 * Falls ein Pfad zu einem Ordner korrekt übergeben wird, dann werden alle
	 * Unterordner als relativer Pfad zurückgegeben.
	 * 
	 * @param pfad
	 *            bis zum Ordner, der die gewünschten Ordner enthält,
	 * @return ArrayList<String> gibt liste von Dateipfaden zurück
	 */
	public static ArrayList<String> ordnerAuslesen(String pfad) {
		File file = new File(pfad);
		ArrayList<String> dateien = new ArrayList<String>();
		if (file.isDirectory()) {
			File[] files = file.listFiles();
			for (int i = 0; i < files.length; i++) {
				try {
					if(files[i].isDirectory())
					dateien.add(files[i].getCanonicalPath());
				} catch (IOException e) {
					System.out.println("Datei konnte nicht gelesen werden");
				}
			}
		}
		return dateien;
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
	 * Either creates a new file or overwrites the existing File with a emtpy
	 * file, and opens the new File for Data input. Data can now be written to
	 * the file via writeLine(). The File must then be closed with closeFile();
	 * 
	 */
	private static void createFile(String pfad) throws IOException,
			FileNotFoundException {
		File file = new File(pfad);
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
			throw new IOException("Cannot create file.");
		}

		// Cancel if File cannot be written to.
		if (!file.canWrite()) {
			throw new IOException("Cannot write to file.");
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
	 * Liest die Datei in eine Hashmap ein.
	 * 
	 * @return Alle Schlüssel - Wert Paare als Hashmap
	 */
	private static HashMap<String, String> readStrings(String pfad)
			throws FileNotFoundException, IOException {

		HashMap<String, String> result = new HashMap<String, String>();
		BufferedReader reader = null;

		try {
			reader = new BufferedReader(new FileReader(pfad));

		} catch (IOException e) {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException ex) {
			}

			throw new FileNotFoundException("File not Found");
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
	 * @param line
	 *            The Line to be written to the File
	 */
	private static void writeLine(String line) throws IOException {
		if (writer == null)
			throw new IOException("Cannot write to file, file not open.");
		writer.println(line);
	}

}