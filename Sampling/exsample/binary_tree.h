// -*- C++ -*-
//
// binary_tree.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2017 Simon Platzer -- simon.plaetzer@desy.de, The Herwig Collaboration
//
// ExSample is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_binary_tree_h_included
#define EXSAMPLE_binary_tree_h_included

#include "utility.h"

namespace exsample {

  /// \brief binary_tree represents a binary tree with the ability to
  /// `cascade' visitor objects down the tree
  template<class Value>
  class binary_tree {

  public:

    ///@name type definitions
    //@{

    /// define the object type
    typedef Value value_type;

    //@}

  public:

    ///@name constructors
    //@{

    /// default constructor
    binary_tree()
      : neighbours_(),
	parent_(), value_(),
	children_()
    { }

    /// construct giving key/cell and parent
    binary_tree(const value_type& thevalue,
		binary_tree * theparent = 0)
      : neighbours_(), parent_(theparent),
	value_(new value_type(thevalue)),
	children_()
    { }

    /// binary_tree has a strict ownership; on copying
    /// binary trees ownership is transferred
    binary_tree(const binary_tree& x)
      : neighbours_(x.neighbours_),
	parent_(x.parent_), value_(),
	children_() {
      assert(x.root());
      binary_tree& nc_x = const_cast<binary_tree&>(x);
      value_.swap(nc_x.value_);
      children_.first.swap(nc_x.children_.first);
      children_.second.swap(nc_x.children_.second);
      nc_x.parent_ = 0;
      nc_x.neighbours_.first = 0;
      nc_x.neighbours_.second = 0;
    }

    /// binary_tree has a strict ownership; on copying
    /// binary trees ownership is transferred
    binary_tree& operator=(const binary_tree& x) {
      if (this == &x)
	return *this;
      assert(x.root());
      binary_tree& nc_x = const_cast<binary_tree&>(x);
      value_.swap(nc_x.value_);
      children_.first.swap(nc_x.children_.first);
      children_.second.swap(nc_x.children_.second);
      neighbours_ = x.neighbours_;
      parent_ = x.parent_;
      nc_x.parent_ = 0;
      nc_x.neighbours_.first = 0;
      nc_x.neighbours_.second = 0;
      return *this;
    }

    //@}

    ///@name standard-conforming leaf iterators
    //@{

  public:

    class const_iterator;

    /// iterator
    class iterator {

    public:

      ///@name type definitions for iterator traits
      //@{

      /// define the iterator category
      typedef std::bidirectional_iterator_tag iterator_category;

      /// define the difference_type
      typedef int difference_type;

      /// define the value type
      typedef Value value_type;

      /// define the reference type
      typedef value_type& reference;

      /// define the pointer type
      typedef value_type * pointer;

      //@}

    public:

      ///@name constructors
      //@{

      /// default constructor
      iterator() : pointee(0), post_end(0), pre_begin(0) { }

      /// constructor taking pointee
      iterator(binary_tree * p, std::size_t end = 0)
	: pointee(p), post_end(end), pre_begin(0) { }

      //@

    public:

      ///@name comparisons
      //@{

      /// comparison
      bool operator==(const iterator& x) const {
	return ((pointee == x.pointee) &&
		(post_end == x.post_end) &&
		(pre_begin == x.pre_begin));
      }

      /// comparison
      bool operator!=(const iterator& x) const { return !(*this == x); }

      //@}

    public:

      ///@name derefrence and indirection
      //@{

      /// dereference
      reference operator*() { return pointee->value(); }

      /// indirection
      pointer operator->() { return &**this; }

      /// return reference to the node
      binary_tree& node() { return *pointee; }

      //@}

      /// return raw pointer to the element pointed to
      binary_tree * get() const { return pointee; }

      ///@name biderectional iterator increment/decrements
      //@{

      /// pre-increment
      iterator& operator++() {
	if (post_end) { ++post_end; return *this; }
	if (pre_begin) { --pre_begin; return *this; }
	if(!(pointee->right_neighbour())) { post_end = 1; return *this; }
	pointee = pointee->right_neighbour();
	return *this;
      }

      /// pre-decrement
      iterator& operator--() {
	if (post_end) { --post_end; return *this; }
	if (pre_begin) { ++pre_begin; return *this; }
	if(!(pointee->left_neighbour())) { pre_begin = 1; return *this; }
	pointee = pointee->left_neighbour();
	return *this;
      }

      /// post-increment
      iterator operator++(int) {
	iterator tmp = *this;
	++(*this);
	return tmp;
      }
      
      /// post-decrement
      iterator operator--(int) {
	iterator tmp = *this;
	--(*this);
	return tmp;
      }
      
      //@}

    private:

      /// friend for conversion
      friend class const_iterator;

      /// the node pointed to
      binary_tree * pointee;

      /// the distance from --end() (if above --end())
      std::size_t post_end;

      /// the distance from begin() (if below begin())
      std::size_t pre_begin;

    };

    /// return begin iterator
    iterator begin() { return iterator(left_most()); }

    /// return end iterator
    iterator end() { return iterator(right_most(),1); }

    /// return global begin iterator
    iterator global_begin() { 
      if (!root())
	return parent().global_begin();
      return iterator(left_most());
    }

    /// return global end iterator
    iterator global_end() {
      if (!root())
	return parent().global_end();
      return iterator(right_most(),1);
    }

    /// const_iterator
    class const_iterator {

    public:

      ///@name type definitions for iterator traits
      //@{

      /// define the iterator category
      typedef std::bidirectional_iterator_tag iterator_category;

      /// define the difference type
      typedef int difference_type;

      /// define the value type
      typedef const Value value_type;

      /// define the reference type
      typedef const value_type& reference;

      /// define the pointer type
      typedef const value_type * pointer;

      //@}

    public:

      ///@name constructors
      //@{

      /// default constructor
      const_iterator() : pointee(0), post_end(0), pre_begin(0) { }

      /// constructor taking pointee
      const_iterator(const binary_tree * p, std::size_t end = 0)
	: pointee(p), post_end(end), pre_begin(0) { }

      /// conversion from iterator
      const_iterator(const iterator& x)
	: pointee(x.pointee), post_end(x.post_end), pre_begin(x.pre_begin) { }

      //@}

    public:

      ///@name comparisons
      //@{

      /// comparison
      bool operator==(const const_iterator& x) const { 
	return ((pointee == x.pointee) &&
		(post_end == x.post_end) &&
		(pre_begin == x.pre_begin));
      }

      /// comparison
      bool operator!=(const const_iterator& x) const { return !(*this == x); }

      //@}

    public:

      ///@name dereference and indirection
      //@{

      /// dereference
      reference operator*() const { return pointee->value(); }

      /// indirection
      pointer operator->() const { return &**this; }

      /// return reference to the node
      const binary_tree& node() const { return *pointee; }

      //@}

      ///@name biderectional iterator increment/decrements
      //@{

      /// pre-increment
      const_iterator& operator++() {
	if (post_end) { ++post_end; return *this; }
	if (pre_begin) { --pre_begin; return *this; }
	if(!(pointee->right_neighbour())) { post_end = 1; return *this; }
	pointee = pointee->right_neighbour();
	return *this;
      }

      /// pre-decrement
      const_iterator& operator--() {
	if (post_end) { --post_end; return *this; }
	if (pre_begin) { ++pre_begin; return *this; }
	if(!(pointee->left_neighbour())) { pre_begin = 1; return *this; }
	pointee = pointee->left_neighbour();
	return *this;
      }

      /// post-increment
      const_iterator operator++(int) {
	const_iterator tmp = *this;
	++(*this);
	return tmp;
      }
      
      /// post-decrement
      const_iterator operator--(int) {
	const_iterator tmp = *this;
	--(*this);
	return tmp;
      }

      //@}

    private:

      /// the node pointed to
      const binary_tree * pointee;

      /// the distance from --end() (if above --end())
      std::size_t post_end;

      /// the distance from begin() (if below begin())
      std::size_t pre_begin;

    };

    /// return begin const_iterator
    const_iterator begin() const { return const_iterator(left_most()); }

    /// return end const_iterator
    const_iterator end() const { return const_iterator(right_most(),1); }

    /// return global begin iterator
    const_iterator global_begin() const { 
      if (!root())
	return parent().global_begin();
      return iterator(left_most());
    }

    /// return global end iterator
    const_iterator global_end() const {
      if (!root())
	return parent().global_end();
      return iterator(right_most(),1);
    }

  private:

    /// set the left neighbour
    void left_neighbour(binary_tree * n) { neighbours_.first = n; }

    /// set the right neighbour
    void right_neighbour(binary_tree * n) { neighbours_.second = n; }

    /// get the left neighbour
    binary_tree * left_neighbour() const { return neighbours_.first; }

    /// get the right neighbour
    binary_tree * right_neighbour() const { return neighbours_.second; }

    /// return the left-most leaf
    binary_tree * left_most() {
      if(leaf()) return this;
      return left_child().left_most();
    }

    /// return the right-most leaf
    binary_tree * right_most() {
      if(leaf()) return this;
      return right_child().right_most();
    }

    /// return the left-most leaf
    const binary_tree * left_most() const {
      if(leaf()) return this;
      return left_child().left_most();
    }

    /// return the right-most leaf
    const binary_tree * right_most() const {
      if(leaf()) return this;
      return right_child().right_most();
    }

    /// the iterator is a good friend
    friend class binary_tree<value_type>::iterator;

    /// the iterator is a good friend
    friend class binary_tree<value_type>::const_iterator;

    /// the left and right neighbours of this node
    std::pair<binary_tree*,binary_tree*> neighbours_;

    //@}

  public:

    /// return true, if this node is empty
    bool empty() const { return root() && leaf() && !value_; }

    /// clear this node
    void clear() {
      neighbours_ = std::make_pair<binary_tree*,binary_tree*>(0,0);
      parent_ = 0;
      value_.reset(0);
      if (!leaf()) {
	left_child().clear();
	right_child().clear();
      }
      children_.first.reset(0);
      children_.second.reset(0);
    }

  public:

    /// split this node
    std::pair<iterator,iterator> split(std::pair<value_type,value_type> children) {

      assert(leaf());

      children_.first.reset(new binary_tree(children.first,this));
      children_.second.reset(new binary_tree(children.second,this));

      children_.first->left_neighbour(neighbours_.first);
      children_.first->right_neighbour(children_.second.get());
      children_.second->left_neighbour(children_.first.get());
      children_.second->right_neighbour(neighbours_.second);

      // adjust original neighbours
      
      if(neighbours_.first) {
	neighbours_.first->right_neighbour(children_.first.get());
      }

      if (neighbours_.second) {
	neighbours_.second->left_neighbour(children_.second.get());
      }

      neighbours_.first = 0; neighbours_.second = 0;

      return std::make_pair(iterator(children_.first.get()),iterator(children_.second.get()));

    }

  public:

    /// select using a selector
    template<class Selector>
    iterator select(const Selector& selector) {

      if(leaf()) {
	bool use = selector.use(value());
	if (use) return iterator(this);
	return global_end();
      }

      std::pair<bool,bool> which(selector.use(value(),left_child().value(),right_child().value()));

      assert(!which.first || !which.second);

      if (!which.first && !which.second) {
	return global_end();
      }

      if (which.first) {
	return left_child().select(selector);
      }
      else {
	return right_child().select(selector);
      }

      return global_end();

    }

    /// generate a hash value for the sub-tree
    /// selected by the given selector object
    template<class Selector, unsigned long bits>
    void subtree_hash(const Selector& selector, bit_container<bits>& bhash) {
      bhash = bit_container<bits>();
      unsigned long pos = 0;
      do_subtree_hash<Selector,bits>(selector,bhash,pos);
    }

    /// accumulate values using a binary function
    /// and accessor object
    template<class Accessor, class BinaryOp>
    typename BinaryOp::result_type accumulate(const Accessor& acc,
					      BinaryOp binary_op) const {

      if (!leaf()) {
	return
	  binary_op(left_child().accumulate(acc,binary_op),
		    right_child().accumulate(acc,binary_op));
      }

      return acc.get(value(),true);

    }

    /// accumulate values only from branches
    /// matching a Selector
    template<class Selector, class Accessor, class BinaryOp>
    typename BinaryOp::result_type accumulate(const Selector& selector,
					      const Accessor& acc,
					      BinaryOp binary_op) const {

      if (!leaf()) {
	std::pair<bool,bool> which(selector.use(value(),left_child().value(),right_child().value()));
	assert(which.first || which.second);
	if (which.first && which.second) {
	  return
	    binary_op(left_child().accumulate(selector,acc,binary_op),
		      right_child().accumulate(selector,acc,binary_op));
	} else if (which.first) {
	  return left_child().accumulate(selector,acc,binary_op);
	} else if (which.second) {
	  return right_child().accumulate(selector,acc,binary_op);
	}
      }

      return acc.get(value(),true);

    }

    /// accumulate values using a binary function
    /// and accessor object, storing intermediate
    /// values in nodes
    template<class Accessor, class BinaryOp>
    typename BinaryOp::result_type tree_accumulate(const Accessor& acc,
						   BinaryOp binary_op) {

      if (!leaf()) {
	acc.set(value()) =
	  binary_op(left_child().tree_accumulate(acc,binary_op),
		    right_child().tree_accumulate(acc,binary_op));
	return acc.get(value(),false);
      }

      acc.set(value()) = acc.get(value(),true);
      return acc.get(value(),true);

    }

    /// accumulate values only from branches
    /// matching a Selector
    template<class Selector, class Accessor, class BinaryOp>
    typename BinaryOp::result_type tree_accumulate(const Selector& selector,
						   const Accessor& acc,
						   BinaryOp binary_op) {

      if (!leaf()) {
	std::pair<bool,bool> which(selector.use(value(),left_child().value(),right_child().value()));
	assert(which.first || which.second);
	if (which.first && which.second) {
	  acc.set(value()) =
	    binary_op(left_child().tree_accumulate(selector,acc,binary_op),
		      right_child().tree_accumulate(selector,acc,binary_op));
	} else if (which.first) {
	  acc.set(value()) = left_child().tree_accumulate(selector,acc,binary_op);
	} else if (which.second) {
	  acc.set(value()) = right_child().tree_accumulate(selector,acc,binary_op);
	}
	return acc.get(value(),false);
      }

      acc.set(value()) = acc.get(value(),true);
      return acc.get(value(),true);

    }

    /// forward propagate a visitor to all children nodes
    template<class Visitor>
    void cascade(Visitor visitor) const {
      if (leaf()) {
	visitor.visit(value()); 
	return;
      } else visitor.visit(value(),left_child().value(),right_child().value());
      left_child().cascade(visitor);
      right_child().cascade(visitor);
    }

    /// succesively split using a generator
    template<class Generator>
    void generate(Generator generator) {
      if (root())
	value_.reset(new value_type(generator.root()));
      if (generator.split()) {
	std::pair<iterator,iterator> ch = split(generator.generate(value()));
	ch.first.node().generate(generator);
	ch.second.node().generate(generator);
      }
    }

  public:

    ///@name Public member access
    //@{

    /// return the value held by this node
    value_type& value() { return *value_; }

    /// return the value held by this node
    const value_type& value() const { return *value_; }

    /// return true, if this is the root node
    bool root() const { return !parent_; }

    /// return true, if this node has got children
    bool leaf() const { return !(children_.first.get() && children_.second.get()); }

    //@}

  public:

    ///@name put and get from streams
    //@{

    /// forward visitor writing out the tree to given ostream
    template<class OStream>
    struct ostream_visitor {

      /// construct from ostream reference
      explicit ostream_visitor(OStream& os) : os_(&os), first_time_(true) {}

      /// visit a leaf node
      void visit(const value_type&) {
	(*os_) << "end_branch";
	ostream_traits<OStream>::separator(*os_);
      }

      /// visit a branching
      void visit(const value_type& parent,
		 const value_type& left, const value_type& right) {
	if (first_time_) {
	  (*os_) << "root_node";
	  ostream_traits<OStream>::separator(*os_);
	  parent.put(*os_);
	  first_time_ = false;
	}

	(*os_) << "left_child";
	ostream_traits<OStream>::separator(*os_);
	left.put(*os_);

	(*os_) << "right_child";
	ostream_traits<OStream>::separator(*os_);
	right.put(*os_);

      }

    private:

      /// pointer to the ostream to write to
      OStream* os_;

      /// whether we are at the or not
      bool first_time_;

    };

    /// generator reading binary tree from istream
    template<class IStream>
    struct istream_generator {

      /// construct from istream reference
      explicit istream_generator(IStream& is)
	: is_(&is), children_(), tag_("") {}

      /// copy constructor
      istream_generator(const istream_generator& x)
	: is_(x.is_), children_(), tag_("") {}

      /// read the root node
      value_type root() {

	*is_ >> tag_;
	assert(tag_ == "root_node");
	value_type rnode;
	rnode.get(*is_);

	return rnode;
      }

      /// read children nodes
      bool split() {

	*is_ >> tag_;

	if (tag_ == "end_branch") {
	  return false;
	}

	assert (tag_ == "left_child");
	children_.first.get(*is_);

	*is_ >> tag_;
	assert(tag_ == "right_child");
	children_.second.get(*is_);

	return true;

      }

      /// return the children generated
      std::pair<value_type,value_type> generate(const value_type&) {
	return children_;
      }

      /// initialize a leaf
      void initialize_leaf(const value_type&) {}

    private:

      /// pointer to the istream used
      IStream* is_;

      /// the children currently handled
      std::pair<value_type,value_type> children_;

      /// temporary storage for tags
      std::string tag_;

    };

    /// put to ostream
    template<class OStream>
    void put(OStream& os) const {

      if (empty()) {
	os << "empty";
	ostream_traits<OStream>::separator(os);
	return;
      } else if (root() && leaf()) {
	os << "root_only";
	ostream_traits<OStream>::separator(os);
	value().put(os);
	return;
      } else {
	os << "non_empty";
	ostream_traits<OStream>::separator(os);
      }

      assert(root());
      cascade(ostream_visitor<OStream>(os));

    }

    /// get from istream
    template<class IStream>
    void get(IStream& is) {

      std::string state;
      is >> state;

      if (state == "empty") {
	return;
      }
      if (state == "root_only") {
	value_.reset(new value_type());
	value().get(is);
	return;
      }

      assert(empty());
      generate(istream_generator<IStream>(is));

    }

    //@}


  private:

    /// calculate hash value
    template<class Selector, unsigned long bits>
    void do_subtree_hash(const Selector& selector,
			 bit_container<bits>& current,
			 unsigned long& position,
			 bool selected = true) const {

      if (!leaf()) {

	std::pair<bool,bool> which(false,false);
	if (selected)
	  which = selector.use(value(),left_child().value(),right_child().value());

	current.bit(position,which.first);
	current.bit(position+1,which.second);

	position += 2;

	left_child().do_subtree_hash(selector,current,position,which.first && selected);
	right_child().do_subtree_hash(selector,current,position,which.second && selected);
      }

    }

  private:

    ///@name private member access
    //@{

    /// return the parent of this node
    binary_tree& parent() { assert(parent_); return *parent_; }

    /// return the parent of this node
    const binary_tree& parent() const { assert(parent_); return *parent_; }

    /// return the left child of this node
    binary_tree& left_child() { assert(children_.first.get()); return *children_.first; }

    /// return the left child of this node
    const binary_tree& left_child() const { assert(children_.first.get()); return *children_.first; }

    /// return the right child of this node
    binary_tree& right_child() { assert(children_.second.get()); return *children_.second; }

    /// return the right child of this node
    const binary_tree& right_child() const { assert(children_.second.get()); return *children_.second; }

    //@}

  private:

    /// the parent of this node
    binary_tree * parent_;

    /// the cell held by this node
    std::unique_ptr<value_type> value_;

    /// the children of this node
    std::pair<std::unique_ptr<binary_tree>,
	      std::unique_ptr<binary_tree> > children_;

  };

}

#endif // EXSAMPLE_binary_tree_h_included
