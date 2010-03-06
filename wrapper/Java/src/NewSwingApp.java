import java.awt.BorderLayout;

import javax.swing.*;

import java.io.*;
import java.lang.reflect.*;

/**
 * This code was edited or generated using CloudGarden's Jigloo
 * SWT/Swing GUI Builder, which is free for non-commercial
 * use. If Jigloo is being used commercially (ie, by a corporation,
 * company or business for any purpose whatever) then you
 * should purchase a license for each developer using Jigloo.
 * Please visit www.cloudgarden.com for details.
 * Use of Jigloo implies acceptance of these licensing terms.
 * A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
 * THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
 * LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
 */
public class NewSwingApp extends javax.swing.JFrame {


	/**
	 * Add custom directory to java.library.path of local JVM
	 * @param s local filesystem directory
	 */
	public static void addDir(String s) throws IOException {
		try {
			// This enables the java.library.path to be modified at runtime
			// From a Sun engineer at http://forums.sun.com/thread.jspa?threadID=707176
			//
			Field field = ClassLoader.class.getDeclaredField("usr_paths");
			field.setAccessible(true);
			String[] paths = (String[])field.get(null);
			for (int i = 0; i < paths.length; i++) {
				if (s.equals(paths[i])) {
					return;
				}
			}
			String[] tmp = new String[paths.length+1];
			System.arraycopy(paths,0,tmp,0,paths.length);
			tmp[paths.length] = s;
			field.set(null,tmp);
			System.setProperty("java.library.path", System.getProperty("java.library.path") + File.pathSeparator + s);
		} catch (IllegalAccessException e) {
			throw new IOException("Failed to get permissions to set library path");
		} catch (NoSuchFieldException e) {
			throw new IOException("Failed to get field handle to set library path");
		}
	}


	/**
	 * Copy library file from jar archive to temp directory of local filesystem
	 * @param name core name of the library 
	 */
	private static String downloadLibrary(String name)
	{
		try
		{
			// Get input stream from jar resource
			// Copy resource to filesystem in a temp folder with a unique name
			File temporaryFile = File.createTempFile("tmpfile", ".file");
			String tmpdir = temporaryFile.getParent();

			InputStream inputStream;
			String path;
			String OS = System.getProperty("os.name").toLowerCase();
			if (OS.indexOf("windows") > -1) {
				inputStream = NewSwingApp.class.getResource("/" + name + ".dll").openStream();  
				path = tmpdir + File.separator + name + ".dll";

			} else {
				if( name.indexOf("lib") < 0){
					name = "lib" + name;
				}
				inputStream = NewSwingApp.class.getResource("/" + name + ".so").openStream();  
				path = tmpdir + File.separator + name + ".so";

			}

			File temporaryDll = new File(path);
			FileOutputStream outputStream = new FileOutputStream(temporaryDll);
			byte[] array = new byte[8192];
			int read = 0;
			while ( (read = inputStream.read(array)) > 0)
				outputStream.write(array, 0, read);
			outputStream.close();  

			// Delete on exit the dll
			temporaryDll.deleteOnExit();  
			temporaryFile.deleteOnExit();  

			System.out.println("saving... " + temporaryDll.getPath());

			// Finally, load the dll
			return new String(temporaryDll.getPath());
		}
		catch(Throwable e)
		{
			e.printStackTrace();
			return new String("");
		}  
	}

	private JMenu jMenu5;
	private JTextField jTextField1;
	private JPanel jPanel1;
	private JMenuItem deleteMenuItem;
	private JSeparator jSeparator1;
	private JMenuItem pasteMenuItem;
	private JMenuItem copyMenuItem;
	private JMenuItem cutMenuItem;
	private JMenu jMenu4;
	private JMenuItem exitMenuItem;
	private JSeparator jSeparator2;
	private JMenuItem closeFileMenuItem;
	private JMenuItem saveAsMenuItem;
	private JMenuItem saveMenuItem;
	private JMenuItem openFileMenuItem;
	private JMenuItem newFileMenuItem;
	private JMenu jMenu3;
	private JMenuBar jMenuBar1;

	static {

		try{ 
			// Find path to temporary directory
			File temporaryDll = File.createTempFile("temp", ".file");
			String tmpdir = temporaryDll.getParent();
			// Add temporary directory to java.lib.path
			addDir(tmpdir);

		} catch(Exception e){
		}

		downloadLibrary("libgslcblas");
		downloadLibrary("libgsl");
		downloadLibrary("example");

		String OS = System.getProperty("os.name").toLowerCase();
		if (OS.indexOf("windows") > -1) {
			System.loadLibrary("libgslcblas");
			System.loadLibrary("libgsl");
		} else {
			System.loadLibrary("gslcblas");
			System.loadLibrary("gsl");
		}

		System.loadLibrary("example");

	}
	
	//public static native int AT_GetNumber();

	/**
	 * Auto-generated main method to display this JFrame
	 */
	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				int res = example.AT_GetNumber();
				NewSwingApp inst = new NewSwingApp(Integer.toString(res));
				inst.setLocationRelativeTo(null);
				inst.setVisible(true);
			}
		});
	}

	public NewSwingApp(String s) {
		super();
		this.initGUI(s);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	private void initGUI(String s) {
		try {
			{
				jPanel1 = new JPanel();
				getContentPane().add(jPanel1, BorderLayout.CENTER);
				{
					jTextField1 = new JTextField();
					jPanel1.add(jTextField1);
					jTextField1.setText(s);
					jTextField1.setPreferredSize(new java.awt.Dimension(155, 83));
				}
			}
			setSize(400, 300);
			{
				jMenuBar1 = new JMenuBar();
				setJMenuBar(jMenuBar1);
				{
					jMenu3 = new JMenu();
					jMenuBar1.add(jMenu3);
					jMenu3.setText("File");
					{
						newFileMenuItem = new JMenuItem();
						jMenu3.add(newFileMenuItem);
						newFileMenuItem.setText("New");
					}
					{
						openFileMenuItem = new JMenuItem();
						jMenu3.add(openFileMenuItem);
						openFileMenuItem.setText("Open");
					}
					{
						saveMenuItem = new JMenuItem();
						jMenu3.add(saveMenuItem);
						saveMenuItem.setText("Save");
					}
					{
						saveAsMenuItem = new JMenuItem();
						jMenu3.add(saveAsMenuItem);
						saveAsMenuItem.setText("Save As ...");
					}
					{
						closeFileMenuItem = new JMenuItem();
						jMenu3.add(closeFileMenuItem);
						closeFileMenuItem.setText("Close");
					}
					{
						jSeparator2 = new JSeparator();
						jMenu3.add(jSeparator2);
					}
					{
						exitMenuItem = new JMenuItem();
						jMenu3.add(exitMenuItem);
						exitMenuItem.setText("Exit");
					}
				}
				{
					jMenu4 = new JMenu();
					jMenuBar1.add(jMenu4);
					jMenu4.setText("Edit");
					{
						cutMenuItem = new JMenuItem();
						jMenu4.add(cutMenuItem);
						cutMenuItem.setText("Cut");
					}
					{
						copyMenuItem = new JMenuItem();
						jMenu4.add(copyMenuItem);

						copyMenuItem.setText(s);
					}
					{
						pasteMenuItem = new JMenuItem();
						jMenu4.add(pasteMenuItem);
						pasteMenuItem.setText("Paste");
					}
					{
						jSeparator1 = new JSeparator();
						jMenu4.add(jSeparator1);
					}
					{
						deleteMenuItem = new JMenuItem();
						jMenu4.add(deleteMenuItem);
						deleteMenuItem.setText("Delete");
					}
				}
				{
					jMenu5 = new JMenu();
					jMenuBar1.add(jMenu5);
					jMenu5.setText("Help");
					{
						final JMenuItem helpMenuItem = new JMenuItem();
						jMenu5.add(helpMenuItem);
						helpMenuItem.setText("Help");
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
