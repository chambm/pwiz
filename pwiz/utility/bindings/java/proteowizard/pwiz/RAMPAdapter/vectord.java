/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.39
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package proteowizard.pwiz.RAMPAdapter;

public class vectord {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected vectord(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(vectord obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if(swigCPtr != 0 && swigCMemOwn) {
      swigCMemOwn = false;
      pwiz_swigbindingsJNI.delete_vectord(swigCPtr);
    }
    swigCPtr = 0;
  }

  public vectord() {
    this(pwiz_swigbindingsJNI.new_vectord__SWIG_0(), true);
  }

  public vectord(long n) {
    this(pwiz_swigbindingsJNI.new_vectord__SWIG_1(n), true);
  }

  public long size() {
    return pwiz_swigbindingsJNI.vectord_size(swigCPtr, this);
  }

  public long capacity() {
    return pwiz_swigbindingsJNI.vectord_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    pwiz_swigbindingsJNI.vectord_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return pwiz_swigbindingsJNI.vectord_isEmpty(swigCPtr, this);
  }

  public void clear() {
    pwiz_swigbindingsJNI.vectord_clear(swigCPtr, this);
  }

  public void add(double x) {
    pwiz_swigbindingsJNI.vectord_add(swigCPtr, this, x);
  }

  public double get(int i) {
    return pwiz_swigbindingsJNI.vectord_get(swigCPtr, this, i);
  }

  public void set(int i, double val) {
    pwiz_swigbindingsJNI.vectord_set(swigCPtr, this, i, val);
  }

}
