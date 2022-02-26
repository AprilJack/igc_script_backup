import java.util.LinkedList;


public class BinaryTree<E extends Comparable<E> & Filterable> {
	private BinaryTreeNode<E> root;

	public void add(E d)
	{
		if(root != null)
			root.add(d);
		else
			root = new BinaryTreeNode<E>(d,null);
	}
	
	
	public void remove(E d)
	{
		if(root != null)
		{
			if(root.getValue().compareTo(d) == 0)
			{

					if(root.hasLeft())
					{
						BinaryTreeNode<E> right = root.getRight();
						root = root.getLeft();
						root.setParent(null);
						BinaryTreeNode<E> most = root.getMost();
						most.setRight(right);
						if(right != null)
						{
							right.setParent(most);
						}
					}
					else if(root.hasRight())
					{
						BinaryTreeNode<E> left = root.getLeft();
						root = root.getRight();
						root.setParent(null);
						BinaryTreeNode<E> least = root.getLeast();
						least.setLeft(left);
						if(left != null)
						{
							left.setParent(least);
						}
					}
					else
						root = null;
			}
			else root.remove(d);
		}
	}

	public void inorder(LinkedList<E> list)
	{
		if(root != null)
			root.inorder(list);
	}
	
	public void postorder(LinkedList<E> list)
	{
		//System.gc();
		if(root != null)
			root.postorder(list);
	}
	
	public void inorder(LinkedList<E> list, E min, E max) {
		if(root != null)
		{
	//		System.gc();
			root.inorder(list,min,max);
		}
		
	}

	public void removeAll() {
		root = null;
		
	}
	/*
	public int filterExons()
	{
		int filtered = 0;
		if(root != null)
		{
			while(root != null)
			{
		//		if(root.getData().filter())
		//			remove(root.getData());
		//		else
		//		{
		///			boolean removed = root.filter(0);
		//			if(!removed) break;
		//		}
		//		filtered++;
			}
		}
		return filtered;
	}
	*/
	public E removeNextInOrder()
	{
		if(root != null)
		{
			if(root.getLeft() != null)
				return root.removeNextInOrder();
			else
			{
				E result = root.getData();
				remove(result);
				return result;
			}
		}
		else
			return null;
	}
	
	public E removeNextPostOrder()
	{
		if(root != null)
		{
			if(root.getRight() != null)
				return root.removeNextPostOrder();
			else
			{
				E result = root.getData();
				remove(result);
				return result;
			}
		}
		else
			return null;
	}
	
	public void print()
	{
		if(root !=null)
			root.print();
	}

}
