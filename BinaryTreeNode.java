
import java.util.*;

public class BinaryTreeNode<E extends Comparable<E> & Filterable> {
	private BinaryTreeNode<E> left, right;	// children; can be null
	private BinaryTreeNode<E> parent;
	private E data;

	public E getData()
	{
		return data;
	}
	
	public void setParent(BinaryTreeNode<E> node)
	{
		parent = node;
	}
	
	
	/**
	 * Constructs leaf node -- left and right are null
	 */
	public BinaryTreeNode(E data, BinaryTreeNode<E> parent) {
		this.data = data; this.left = null; this.right = null;
		this.parent =parent;
	}
	
	public BinaryTreeNode()
	{
		this.data = null; this.left = null; this.right = null;
	}

	public void add(E d)
	{
		if(this.data == null)
			this.data = d;
		else
		{
			if(d.compareTo(data) > 0)
			{
				if(right  == null)
					right = new BinaryTreeNode<E>(d,this);
				else
					right.add(d);
			}
			else
			{
				if(left == null)
					left = new BinaryTreeNode<E>(d,this);
				else
					left.add(d);
			}
		}
	}
	
	/**
	 * Will only be called if the removed node is not the root!
	 * @param d
	 */
	public void remove(E d)
	{
		if(d.compareTo(data) == 0)
		{
			if(data.compareTo(parent.data) > 0)
			{
				if(left != null)
				{
					parent.right = left;
					left.parent = parent;
					BinaryTreeNode<E> most = left.getMost();
					most.right =right;
					if(right != null)
						right.parent = most;
				}
				else
				{
					//right must exist
					parent.right = right;
					if(right != null)
						right.parent = parent;
				}
			}
			else
			{
				//its to the left
				if(right != null )
				{
					parent.left = right;
					right.parent = parent;
					if(left != null)
					{
						BinaryTreeNode<E> least = right.getLeast();
						least.left =left;
						left.parent = least;
					}
				}
				else
				{
					//left must exist
					parent.left = left;
					if(left != null)
						left.parent = parent;
				}
			}
		}
		else if(d.compareTo(data) < 0)
		{
			if(left != null && !left.equals(this))
				left.remove(d);
		}
		else if(right != null && !right.equals(this))
			right.remove(d);
	}
	
	public BinaryTreeNode<E> getMost()
	{
		if(right == null || right.equals(this))
			return this;
		else
			return right.getMost();
	}
	
	public BinaryTreeNode<E> getLeast()
	{
		if(left == null || left.equals(this))
			return this;
		else
			return left.getLeast();
	}
	

	/**
	 * Does it have a left child?
	 */
	public boolean hasLeft() {
		return left != null;
	}

	/**
	 * Does it have a right child?
	 */
	public boolean hasRight() {
		return right != null;
	}

	/**
	 * @return its left child 
	 */
	public BinaryTreeNode<E> getLeft() {
		return left;
	}

	/**
	 * @return its right child 
	 */
	public BinaryTreeNode<E> getRight() {
		return right;
	}

	/**
	 * Sets its left child to be newLeft 
	 * @param newLeft - the new left child
	 */
	public void setLeft(BinaryTreeNode<E> newLeft) {
		left = newLeft;
	}

	/**
	 * Sets its right child to be newRight 
	 * @param newRight - the new right child
	 */
	public void setRight(BinaryTreeNode<E> newRight) {
		right = newRight;
	}

	/**
	 * @return its data value 
	 */
	public E getValue() {
		return data;
	}

	/**
	 * Sets its data value 
	 * @param newValue - the new data value
	 */
	public void setValue(E newValue) {
		data = newValue;
	}


	/**
	 * Same structure and data?
	 */
	public boolean equals(Object other) {
		if(!(other instanceof BinaryTreeNode<?>))
			return false;
		BinaryTreeNode<E> t2 = (BinaryTreeNode<E>) other;
		if (hasLeft() != t2.hasLeft() || hasRight() != t2.hasRight()) return false;
		if (!data.equals(t2.data)) return false;
		if (hasLeft() && !left.equals(t2.left)) return false;
		if (hasRight() && !right.equals(t2.right)) return false;
		return true;
	}



	/**
	 * Returns a string representation of the tree
	 */
	public String toString() {
		return toStringHelper("");
	}

	/**
	 * Recursively constructs a String representation of the tree from this node, 
	 * starting with the given indentation and indenting further going down the tree
	 */
	private String toStringHelper(String indent) {
		String ret = "";
		if (hasRight()) 
			ret += right.toStringHelper(indent+"  ");
		ret += indent + data + "\n";
		if (hasLeft()) 
			ret += left.toStringHelper(indent+"  ");
		return ret;
	}

	/** 
	 * Creates a list storing the the data values in the subtree of a node,
	 * ordered according to the preorder traversal of the subtree. 
	 * @param dataList - the list to be returned.
	 */
	public void preorder(List<E> dataList) {
		dataList.add(data);
		if (this.hasLeft())
			left.preorder(dataList);	// recurse on left child
		if (this.hasRight())
			right.preorder(dataList);	// recurse on right child
	}

	/** 
	 * Creates a list storing the the data values in the subtree of a node,
	 * ordered according to the inorder traversal of the subtree. 
	 * @param dataList - the list to be returned.
	 */
	public void inorder(List<E> dataList) {
		try{
			if (this.hasLeft() && !left.data.equals(data))
				left.inorder(dataList);	// recurse on left child
			dataList.add(data);
			if (this.hasRight() && !right.data.equals(data))
				right.inorder(dataList);	// recurse on right child
		}catch (StackOverflowError er)
		{
			er.printStackTrace();
		}
	}
	
	public void inorder(List<E> dataList, E min, E max) {
		if (this.hasLeft() && min.compareTo(data) <= -1 && !left.data.equals(data))
			left.inorder(dataList,min,max);
		if(data.compareTo(min) >= -1 && data.compareTo(max) <= 1 ) 
			dataList.add(data);
		if(this.hasRight() && max.compareTo(data) >= 1 && !right.data.equals(data))
			right.inorder(dataList,min,max);		
	}
	

	/** 
	 * Creates a list storing the the data values in the subtree of a node,
	 * ordered according to the postorder traversal of the subtree. 
	 * @param dataList - the list to be returned.
	 */
	public void postorder(List<E> dataList) {

		if (this.hasRight() && !right.data.equals(data))
			right.postorder(dataList);	// recurse on right child
		dataList.add(data);
		if (this.hasLeft() && !left.data.equals(data))
			left.postorder(dataList);	// recurse on left child
	}
	
	public E removeNextInOrder()
	{
		if(left == null)
		{
			E result = data;
			remove(data);
			return result;
		}
		else
			return left.removeNextInOrder();
	}
	
	public E removeNextPostOrder()
	{
		if(right == null)
		{
			E result = data;
			remove(data);
			return result;
		}
		else
			return right.removeNextPostOrder();
	}
	/*
	public boolean filter(int depth) {
	//	if(depth > 1000)
	//		return false;
		if(data.filter())
		{
			remove(data);
			return true;
		}
		else //data is ok
		{
			if(left != null&& !left.equals(this))
			{
				if(left.filter(depth+1))
					return true;
			}
			if(right != null && !right.equals(this))
			{
				if(right.filter(depth+1))
					return true;
			}
			return false;
		}
	}
*/
	public void print() {
		if(left != null) left.print();
		System.out.println(data.toString());
		if(right != null) right.print();
		
	}


	
}