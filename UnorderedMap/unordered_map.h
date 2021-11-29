#include<iostream>
#include<memory>
#include<vector>

template <typename T, typename Allocator = std::allocator<T>>
class List {
 public:
  struct Node {
    Node* next;
    Node* prev;
    T x;
    explicit Node() : next(nullptr), prev(nullptr), x(T()) {}
    explicit Node(const T& x) : next(nullptr), prev(nullptr), x(x) {}
    Node(T&& y) : next(nullptr), prev(nullptr), x(std::move(y)) {}
    Node(Node&& y) : next(y.prev), prev(y.next), x(std::move(y.x)) {}
  };
  size_t length;
  Node* term;
  Allocator alloc;
  using AllocTraits = std::allocator_traits<Allocator>;
  using NodeAllocType = typename AllocTraits::template rebind_alloc<Node>;
  NodeAllocType NodeAlloc;
  using NodeAllocTraits = std::allocator_traits<NodeAllocType>;

  template<bool is_const>
  class TemplateIterator;

  Node* link(Node* node, const T& x) {
    Node* newNode = NodeAllocTraits::allocate(NodeAlloc, 1);
    NodeAllocTraits::construct(NodeAlloc, newNode, x);
    newNode->next = node;
    node->prev->next = newNode;
    newNode->prev = node->prev;
    node->prev = newNode;
    ++length;
    return newNode;
  }

  Node* unlink(Node* node) {
    node->next->prev = node->prev;
    node->prev->next = node->next;
    Node* next = node->next;
    --length;
    NodeAllocTraits::destroy(NodeAlloc, node);
    NodeAllocTraits::deallocate(NodeAlloc, node, 1);
    return next;
  }
 public:
  using iterator = TemplateIterator<false>;
  using const_iterator = TemplateIterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  explicit List(const Allocator& alloc = Allocator()) : alloc(alloc) {
    length = 0;
    NodeAllocType newNodeAlloc;
    NodeAlloc = newNodeAlloc;
    term = NodeAllocTraits::allocate(NodeAlloc, 1);
    term->next = term;
    term->prev = term;
  }

  List(size_t count, const T& x, const Allocator& alloc = Allocator()) : alloc(alloc) {
    length = count;
    NodeAllocType newNodeAlloc;
    NodeAlloc = newNodeAlloc;
    term = NodeAllocTraits::allocate(NodeAlloc, 1);
    term->next = term;
    term->prev = term;
    for (size_t i = 0; i < count; ++i) {
      push_back(x);
    }
  }

  List(size_t count, const Allocator& alloc = Allocator()) : alloc(alloc) {
    length = count;
    NodeAllocType newNodeAlloc;
    NodeAlloc = newNodeAlloc;
    term = NodeAllocTraits::allocate(NodeAlloc, 1);
    term->next = term;
    term->prev = term;
    Node* cur = term;
    for (size_t i = 0; i < count; ++i) {
      Node* new_element = NodeAllocTraits::allocate(NodeAlloc, 1);
      NodeAllocTraits::construct(NodeAlloc, new_element);
      new_element->prev = cur;
      cur->next = new_element;
      new_element->next = term;
      term->prev = new_element;
      cur = new_element;
    }
  }

  Allocator get_allocator() const {
    return alloc;
  }

  List(const List& x) {
    alloc = AllocTraits::select_on_container_copy_construction(x.alloc);
    NodeAlloc = NodeAllocTraits::select_on_container_copy_construction(x.NodeAlloc);
    Node* curOther = x.term->next;
    term = NodeAllocTraits::allocate(NodeAlloc, 1);
    term->next = term;
    term->prev = term;
    length = 0;
    for (size_t i = 0; i < x.length; ++i) {
      push_back(curOther->x);
      curOther = curOther->next;
    }
  }

  List(List&& x) {
    length = x.length;
    alloc = std::move(x.alloc);
    NodeAlloc = std::move(x.NodeAlloc);
    Node* nullterm = NodeAllocTraits::allocate(NodeAlloc, 1);
    nullterm->next = nullterm;
    nullterm->prev = nullterm;
    term = x.term;
    x.term = nullterm;
    x.length = 0;
  }

  ~List() {
    Node* cur = term->next;
    while (cur != term) {
      Node* next = cur->next;
      NodeAllocTraits::destroy(NodeAlloc, cur);
      NodeAllocTraits::deallocate(NodeAlloc, cur, 1);
      cur = next;
    }
    NodeAllocTraits::deallocate(NodeAlloc, cur, 1);
  }

  List& operator=(const List& x) {
    if (this == &x)
      return *this;
    if (AllocTraits::propagate_on_container_copy_assignment::value && alloc != x.alloc) {
      alloc = x.alloc;
      NodeAlloc = x.NodeAlloc;
    }
    List copy = x;
    std::swap(copy.length, this->length);
    std::swap(copy.term, this->term);
    return *this;
  }

  List& operator=(List&& x) {
    if (this == &x)
      return *this;
    if (AllocTraits::propagate_on_container_copy_assignment::value && alloc != x.alloc) {
      alloc = std::move(x.alloc);
      NodeAlloc = std::move(x.NodeAlloc);
    }
    List copy = std::move(x);
    std::swap(copy.length, this->length);
    std::swap(copy.term, this->term);
    return *this;
  }

  iterator insert(const_iterator it, const T& x) {
    Node* iter = const_cast<Node*>(it.cur);
    return iterator(link(iter, x));
  }

  iterator erase(const_iterator it) {
    Node* iter = const_cast<Node*>(it.cur);
    return iterator(unlink(iter));
  }

  size_t size() const {
    return length;
  }

  void push_back(const T& x) {
    link(term, x);
  }

  void push_front(const T& x) {
    link(term->next, x);
  }

  void pop_back() {
    unlink(term->prev);
  }

  void pop_front() {
    unlink(term->next);
  }

  iterator begin() const {
    return iterator(term->next);
  }

  iterator end() const {
    return iterator(term);
  }

  const_iterator cbegin() const {
    return const_iterator(term->next);
  }

  const_iterator cend() const {
    return const_iterator(term);
  }

  reverse_iterator rbegin() const {
    return reverse_iterator(term);
  }

  reverse_iterator rend() const {
    return reverse_iterator(term->next);
  }

  const_reverse_iterator crbegin() const {
    return const_reverse_iterator(term);
  }

  const_reverse_iterator crend() const {
    return const_reverse_iterator(term->next);
  }

  void swap(List& other) {
    std::swap(term, other.term);
    std::swap(length, other.length);
    std::swap(alloc, other.alloc);
    std::swap(NodeAlloc, other.NodeAlloc);
  }
};

template<typename T, typename Allocator>
template<bool is_const>
class List<T, Allocator>::TemplateIterator {
 public:
  using type = Node*;
  type cur;
 public:
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = std::conditional_t<is_const, const T, T>;
  using reference = std::conditional_t<is_const, const T&, T&>;
  using pointer = std::conditional_t<is_const, const T*, T*>;
  using difference_type = std::ptrdiff_t;

  TemplateIterator(Node* cur) : cur(cur) {}
  TemplateIterator(const TemplateIterator& cur) : cur(cur.cur) {}

  operator TemplateIterator<1>() {
    return TemplateIterator<1>(cur);
  }
  bool operator==(const TemplateIterator& x) const {
    return cur == x.cur;
  }
  bool operator!=(const TemplateIterator& x) const {
    return cur != x.cur;
  }
  std::conditional_t<is_const, const T&, T&> operator*() const {
    return cur->x;
  }
  std::conditional_t<is_const, const T*, T*> operator->() const {
    return &(cur->x);
  }
  TemplateIterator& operator++() {
    cur = cur->next;
    return *this;
  }
  TemplateIterator& operator--() {
    cur = cur->prev;
    return *this;
  }
  TemplateIterator operator++(int) {
    TemplateIterator copy = *this;
    cur = cur->next;
    return copy;
  }
  TemplateIterator operator--(int) {
    TemplateIterator copy = *this;
    cur = cur->prev;
    return copy;
  }
};

template<typename Key, typename Value, typename Hash = std::hash<Key>,
    typename Equal = std::equal_to<Key>, typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
 private:
  Hash hash_function;
  Equal equal_function;
  template<bool is_const>
  class TemplateIterator;

  using iterator = TemplateIterator<0>;
  using const_iterator = TemplateIterator<1>;
 public:
  using NodeType = std::pair<const Key, Value>;
  static const size_t minimum_size = 5;
 private:
  using ListAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<NodeType*>;
  std::vector<typename List<NodeType*, ListAlloc>::iterator> buckets;
  std::vector<typename List<NodeType*, ListAlloc>::iterator> ends_buckets;
  List<NodeType*, ListAlloc> list_;
  Alloc alloc;

  float cur_max_load_factor = 1;

  size_t get_bucket_idx(const Key& x) const {
    return hash_function(x) % buckets.size();
  }

  void try_to_expand() {
    if (load_factor() > cur_max_load_factor) {
      rehash(3 * size() + 5);
    }
  }
 public:
  void max_load_factor(float new_load_factor) {
    cur_max_load_factor = new_load_factor;
  }

  float max_load_factor() const {
    return cur_max_load_factor;
  }

  double load_factor() const {
    return static_cast<double>(list_.size()) / static_cast<double>(buckets.size());
  }

  size_t size() const {
    return list_.size();
  }

  UnorderedMap() {
    buckets.resize(minimum_size, list_.end());
  }

  UnorderedMap(const UnorderedMap& x) {
    alloc = std::allocator_traits<Alloc>::select_on_container_copy_construction(x.alloc);
    cur_max_load_factor = x.cur_max_load_factor;
    equal_function = x.equal_function;
    hash_function = x.hash_function;
    buckets.resize(minimum_size, list_.end());
    ends_buckets.resize(minimum_size, list_.end());
    for (auto& element : x) {
      insert(element);
    }
  }

  UnorderedMap(UnorderedMap&& x) {
    list_ = std::move(x.list_);
    buckets = std::move(x.buckets);
    ends_buckets = std::move(x.ends_buckets);
    equal_function = std::move(x.equal_function);
    hash_function = std::move(x.hash_function);
    alloc = std::move(std::allocator_traits<Alloc>::select_on_container_copy_construction(x.alloc));
    cur_max_load_factor = std::move(x.cur_max_load_factor);
  }

  UnorderedMap& operator=(const UnorderedMap& x) {
    if (this == &x) {
      return *this;
    }
    if (std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value && alloc != x.alloc) {
      alloc = x.alloc;
    }
    hash_function = x.hash_function;
    equal_function = x.equal_function;
    buckets.resize(minimum_size, list_.end());
    ends_buckets.resize(minimum_size, list_.end());
    cur_max_load_factor = x.cur_max_load_factor;
    for (auto element : x) {
      insert(element);
    }
    return *this;
  }

  UnorderedMap& operator=(UnorderedMap&& x) {
    if (this == &x) {
      return *this;
    }
    if (std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value && alloc != x.alloc) {
      alloc = std::move(x.alloc);
    }
    buckets = std::move(x.buckets);
    ends_buckets = std::move(x.ends_buckets);
    list_ = std::move(x.list_);
    equal_function = std::move(x.equal_function);
    hash_function = std::move(x.hash_function);
    cur_max_load_factor = std::move(x.cur_max_load_factor);
    return *this;
  }

  ~UnorderedMap() {
    for (auto element : list_) {
      delete element;
    }
  }

  iterator begin() {
    return iterator(list_.begin());
  }

  const_iterator cbegin() const {
    return const_iterator(list_.begin());
  }

  const_iterator begin() const {
    return const_iterator(list_.begin());
  }

  iterator end() {
    return iterator(list_.end());
  }

  const_iterator cend() const {
    return const_iterator(list_.end());
  }

  const_iterator end() const {
    return const_iterator(list_.end());
  }

  iterator find(const Key& x) const {
    size_t index = get_bucket_idx(x);
    iterator it = buckets[index];
    for (; it != list_.end(); ++it) {
      if (equal_function(it->first, x)) {
        return it;
      }
      if (it == ends_buckets[index])
        break;
    }
    return list_.end();
  }

  void reserve(size_t newsize) {
    if (newsize * buckets.size() > size()) {
      rehash(newsize);
    }
  }

  void rehash(size_t newsize) {
    buckets.clear();
    ends_buckets.clear();
    auto cp = std::move(list_);
    buckets.resize(newsize, list_.end());
    ends_buckets.resize(newsize, list_.end());
    for (auto it = cp.begin(); it != cp.end(); ++it) {
      size_t index = get_bucket_idx((*it)->first);
      auto new_element = buckets[index];
      if (buckets[index] == list_.end()) {
        ends_buckets[index] = new_element;
      }
      new_element = list_.insert(new_element, *it);
    }
  }

  Value& at(const Key& x) {
    iterator it = find(x);
    if (it != end()) {
      return it->second;
    }
    throw std::out_of_range("Error, this element has not found");
  }

  const Value& at(const Key& x) const {
    iterator it = find(x);
    if (it != end()) {
      return it->second;
    }
    throw std::out_of_range("Error, this element has not found");
  }

  Value& operator[](const Key& x) {
    try {
      return at(x);
    } catch (...) {
      return insert(std::make_pair(x, Value())).first->second;
    }
  }

  template<typename ...Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    NodeType* new_element = std::allocator_traits<Alloc>::allocate(alloc, 1);
    std::allocator_traits<Alloc>::construct(alloc, new_element, std::forward<Args>(args)...);

    iterator it = find(new_element->first);
    if (it != end()) {
      return std::make_pair(it, 0);
    }
    try_to_expand();

    size_t index = get_bucket_idx(new_element->first);
    auto first_in_cell = buckets[index];
    auto res = list_.insert(first_in_cell, new_element);

    buckets[index] = res;
    if (first_in_cell == list_.end()) {
      ends_buckets[index] = res;
    }

    return std::make_pair(iterator(res), 1);
  }

  std::pair<iterator, bool> insert(const NodeType& x) {
    return emplace(x);
  }

  template<typename G>
  std::pair<iterator, bool> insert(G&& x) {
    return emplace(std::forward<G>(x));
  }

  template<typename InputIterator>
  void insert(InputIterator b, InputIterator e) {
    for (auto i = b; i != e; ++i) {
      insert(*i);
    }
  }

  void erase(const_iterator b, const_iterator e) {
    while (b != e) {
      auto next_element = b;
      ++next_element;
      erase(b);
      b = next_element;
    }
  }

  void erase(const_iterator x) {
    size_t ind = get_bucket_idx(x->first);
    bool last = (x == static_cast<const_iterator>(ends_buckets[ind]));
    bool first = (x == static_cast<const_iterator>(buckets[ind]));
    auto prev = x.value_->prev;
    auto del = list_.erase(x.value_);
    if (!last && !first) {
      return;
    }
    if (first && last) {
      buckets[ind] = list_.end();
      ends_buckets[ind] = list_.end();
      return;
    }
    if (first) {
      buckets[ind] = del;
      return;
    }
    if (last) {
      ends_buckets[ind] = prev;
    }
  }
};

template<typename Key, typename Value, typename Hash,
    typename Equal, typename Alloc>
template<bool is_const>
class UnorderedMap<Key, Value, Hash, Equal, Alloc>::TemplateIterator {
 private:
  std::conditional_t<is_const, typename List<NodeType*, ListAlloc>::const_iterator, typename List<NodeType*, ListAlloc>::iterator> value_;
  friend class UnorderedMap<Key, Value, Hash, Equal, Alloc>;
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename std::conditional_t<is_const, const NodeType, NodeType>;
  using reference = typename std::conditional_t<is_const, const NodeType&, NodeType&>;
  using pointer = typename std::conditional_t<is_const, const NodeType*, NodeType*>;
  using difference_type = std::ptrdiff_t;

  TemplateIterator(const std::conditional_t<is_const, typename List<NodeType*, ListAlloc>::const_iterator, typename List<NodeType*, ListAlloc>::iterator>& x) : value_(x) {}

  operator TemplateIterator<true>() {
    return TemplateIterator<true>(value_);
  }

  bool operator==(const TemplateIterator& x) const {
    return value_ == x.value_;
  }

  bool operator!=(const TemplateIterator& x) const {
    return value_ != x.value_;
  }

  TemplateIterator& operator++() {
    ++value_;
    return *this;
  }

  TemplateIterator operator++(int) {
    auto copy = *this;
    ++value_;
    return copy;
  }

  reference operator*() const {
    return **value_;
  }

  pointer operator->() const {
    return *value_;
  }
};
