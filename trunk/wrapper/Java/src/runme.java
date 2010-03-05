// runme.java

public class runme {
  static {
    System.loadLibrary("example");
  }

  public static void main(String argv[]) {
    System.out.println(example.fact(4));
    example.AT_PrintName();
    System.out.println(example.AT_GetNumber());
  }
}

