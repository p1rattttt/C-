#include <iostream>
#include <type_traits>
#include <memory>

struct BaseDeleter {
  virtual void operator()(void*) = 0;
  virtual ~BaseDeleter() = default;
};

template<typename T, typename U>
struct DeleterWithAllocator : BaseDeleter {
  U allocator;

  DeleterWithAllocator(U allocator) : allocator(allocator) {}

  void operator()(void* ptr) override {
    using right_allocator_type = typename std::allocator_traits<U>::template rebind_alloc<T>;
    right_allocator_type right_allocator = allocator;
    using right_traits = std::allocator_traits<right_allocator_type>;
    right_traits::destroy(right_allocator, reinterpret_cast<T*>(ptr));
  }
};

template<typename T, typename U = std::default_delete<T>>
struct OwnDeleter : BaseDeleter {
  U deleter;

  OwnDeleter(U deleter) : deleter(deleter) {}

  void operator()(void* ptr) override {
    deleter(reinterpret_cast<T*>(ptr));
  }
};

struct BaseAllocator {
  virtual void dealloc(void* ptr) = 0;
  virtual bool TAllocated() = 0;
  virtual ~BaseAllocator() = default;
};

template<typename T, typename U = std::allocator<char>>
struct AllocatorWithNotCstyle : BaseAllocator {
  U allocator;

  AllocatorWithNotCstyle(U allocator) : allocator(allocator) {}

  void dealloc(void* ptr) override {
    using Right_allocator_type = typename std::allocator_traits<U>::template rebind_alloc<char>;
    Right_allocator_type right_allocator = allocator;
    using right_traits = std::allocator_traits<Right_allocator_type>;
    right_traits::deallocate(right_allocator,
                             reinterpret_cast<char*>(ptr),
                             sizeof(T) + sizeof(BaseDeleter) + sizeof(BaseAllocator) + 2 * sizeof(size_t));
  }

  bool TAllocated() override {
    return true;
  }
};

template<typename T, typename U = std::allocator<char>>
struct AllocatorWithCStyle : BaseAllocator {
  U allocator;

  AllocatorWithCStyle(U allocator) : allocator(allocator) {}

  void dealloc(void* ptr) override {
    using Right_allocator_type = typename std::allocator_traits<U>::template rebind_alloc<char>;
    Right_allocator_type right_allocator = allocator;
    using right_traits = std::allocator_traits<Right_allocator_type>;
    right_traits::deallocate(right_allocator,
                             reinterpret_cast<char*>(ptr),
                             sizeof(BaseDeleter) + sizeof(BaseAllocator) + 2 * sizeof(size_t));
  }

  bool TAllocated() override {
    return false;
  }
};

template<typename T>
class WeakPtr;

template<typename T>
class SharedPtr {
 private:
  template<typename U, typename... Args>
  friend SharedPtr<U> makeShared(Args&& ...);

  template<typename U, typename Allocator, typename... Args>
  friend SharedPtr<U> allocateShared(const Allocator&, Args&& ...);

  template<typename U>
  friend
  class SharedPtr;

  template<class U>
  friend
  class WeakPtr;
// Something like ControlBlock
  T* ptr_ = nullptr;
  BaseDeleter* deleter_ = nullptr;
  BaseAllocator* allocator_ = nullptr;
  size_t* shared_count_ = nullptr;
  size_t* weak_count_ = nullptr;
//
  template<typename U>
  void AllocateBlockWithoutPtr(U allocator) {
    char* ControlBlock =
        std::allocator_traits<U>::allocate(allocator, sizeof(BaseDeleter) + sizeof(BaseAllocator) + 2 * sizeof(size_t));
    deleter_ = reinterpret_cast<BaseDeleter*>(ControlBlock);
    allocator_ = reinterpret_cast<BaseAllocator*>(ControlBlock + sizeof(BaseDeleter));
    shared_count_ = reinterpret_cast<size_t*>(ControlBlock + sizeof(BaseDeleter) + sizeof(BaseAllocator));
    weak_count_ =
        reinterpret_cast<size_t*>(ControlBlock + sizeof(BaseDeleter) + sizeof(BaseAllocator) + sizeof(size_t));
  }

  void destroy() {
    ptr_ = nullptr;
    deleter_ = nullptr;
    allocator_ = nullptr;
    shared_count_ = nullptr;
    weak_count_ = nullptr;
  }

  template<class U>
  void assigne(U&& other) {
    ptr_ = other.ptr_;
    deleter_ = other.deleter_;
    allocator_ = other.allocator_;
    shared_count_ = other.shared_count_;
    weak_count_ = other.weak_count_;
  }

  SharedPtr(const WeakPtr<T>& x) {
    equal(x);
    if (shared_count_)
      ++*shared_count_;
  }

  SharedPtr(char* ControlBlock) { // Constructor for allocateShared
    ptr_ = reinterpret_cast<T*>(ControlBlock);
    deleter_ = reinterpret_cast<BaseDeleter*>(sizeof(T) + ControlBlock);
    allocator_ = reinterpret_cast<BaseAllocator*>(sizeof(T) + ControlBlock + sizeof(BaseDeleter));
    shared_count_ = reinterpret_cast<size_t*>(sizeof(T) + ControlBlock + sizeof(BaseDeleter) + sizeof(BaseAllocator));
    weak_count_ = reinterpret_cast<size_t*>(sizeof(T) + ControlBlock + sizeof(BaseDeleter) + sizeof(BaseAllocator)
        + sizeof(size_t));
  }
 public:
  SharedPtr() = default;

  T* get() {
    return ptr_;
  }

  const T* get() const {
    return ptr_;
  }

  size_t use_count() const {
    return *shared_count_;
  }

  T* operator->() const {
    return ptr_;
  }

  T& operator*() const {
    return *ptr_;
  }

  void reset() {
    SharedPtr<T>().swap(*this);
  }

  template<class U>
  void reset(U* ptr) {
    SharedPtr<T>(ptr).swap(*this);
  }

  template<typename U, typename Deleter>
  void reset(U* ptr, Deleter deleter) {
    SharedPtr<T>(ptr, deleter).swap(*this);
  }

  template<typename U, typename Deleter = std::default_delete<T>, typename Allocator = std::allocator<U>>
  void reset(U* ptr, Deleter deleter, Allocator allocator) {
    SharedPtr<T>(ptr, deleter, allocator).swap(*this);
  }

  template<class U, class Deleter = std::default_delete<T>, class Allocator = std::allocator<U>>
  SharedPtr(U* ptr, Deleter deleter = Deleter(), Allocator allocator = Allocator())
      : ptr_(ptr) { //Constructor for C-Style Pointer
    using AllocateAlloc = typename std::allocator_traits<Allocator>::template rebind_alloc<char>;
    AllocateAlloc allocate_alloc = allocator;
    AllocateBlockWithoutPtr(allocate_alloc);
    *shared_count_ = 1;
    *weak_count_ = 0;
    new(deleter_) OwnDeleter<T, Deleter>(deleter);
    new(allocator_) AllocatorWithCStyle<T, Allocator>(allocator);
  }

  SharedPtr(const SharedPtr& other) {
    equal(other);
    if (shared_count_)
      ++*shared_count_;
  }

  template<class U>
  SharedPtr(const SharedPtr<U>& other) {
    assigne(other);
    if (shared_count_)
      ++*shared_count_;
  }

  SharedPtr(SharedPtr&& other) {
    assigne(std::move(other));
    other.destroy();
  }

  template<class U>
  SharedPtr(SharedPtr<U>&& other) {
    assigne(std::move(other));
    other.destroy();
  }

  SharedPtr& operator=(const SharedPtr& other) {
    SharedPtr copy = other;
    swap(copy);
    return *this;
  }

  template<class U>
  SharedPtr& operator=(const SharedPtr<U>& other) {
    SharedPtr copy = other;
    swap(copy);
    return *this;
  }

  SharedPtr& operator=(SharedPtr&& other) {
    SharedPtr copy = std::move(other);
    swap(copy);
    return *this;
  }

  template<class U>
  SharedPtr& operator=(SharedPtr<U>&& other) {
    SharedPtr copy = std::move(other);
    swap(copy);
    return *this;
  }

  void swap(SharedPtr& x) {
    std::swap(ptr_, x.ptr_);
    std::swap(deleter_, x.deleter_);
    std::swap(allocator_, x.allocator_);
    std::swap(shared_count_, x.shared_count_);
    std::swap(weak_count_, x.weak_count_);
  }

  ~SharedPtr() {
    if (shared_count_ == nullptr) return;
    --*shared_count_;
    if (!(*shared_count_)) {
      deleter_->operator()(ptr_);
      if (!(*weak_count_)) {
        if (allocator_->TAllocated()) { // if ptr was not C-style
          allocator_->dealloc(ptr_);
        } else {
          allocator_->dealloc(deleter_);
        }
      }
    }
  }
};

template<typename T, typename Allocator, typename... Args>
SharedPtr<T> allocateShared(const Allocator& alloc, Args&& ... args) {
  using charAllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<char>;
  using TAllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<T>;
  charAllocatorType charAllocator = alloc;
  TAllocatorType TAllocator = alloc;
  char* ControlBlock = std::allocator_traits<charAllocatorType>::allocate(charAllocator,
                                                                          sizeof(T) + sizeof(BaseDeleter)
                                                                              + sizeof(BaseAllocator)
                                                                              + 2 * sizeof(size_t));
  std::allocator_traits<TAllocatorType>::construct(TAllocator,
                                                   reinterpret_cast<T*>(ControlBlock),
                                                   std::forward<Args>(args)...);
  auto deleter_ = reinterpret_cast<BaseDeleter*>(ControlBlock + sizeof(T));
  auto allocator_ = reinterpret_cast<BaseAllocator*>(ControlBlock + sizeof(T) + sizeof(BaseDeleter));
  auto
      shared_count_ = reinterpret_cast<size_t*>(ControlBlock + sizeof(T) + sizeof(BaseDeleter) + sizeof(BaseAllocator));
  auto weak_count_ = reinterpret_cast<size_t*>(ControlBlock + sizeof(T) + sizeof(BaseDeleter) + sizeof(BaseAllocator)
      + sizeof(size_t));
  *shared_count_ = 1;
  *weak_count_ = 0;
  new(deleter_) DeleterWithAllocator<T, Allocator>(TAllocator);
  new(allocator_) AllocatorWithNotCstyle<T, Allocator>(TAllocator);
  return SharedPtr<T>(ControlBlock);
}

template<typename T, typename... Args>
SharedPtr<T> makeShared(Args&& ... args) {
  return allocateShared<T>(std::allocator<T>(), std::forward<Args>(args)...);
}

template<typename T>
class WeakPtr {
 private:
  template<class U>
  friend
  class WeakPtr;

  template<class U>
  friend
  class SharedPtr;

  //ControlBlock
  T* ptr_ = nullptr;
  BaseDeleter* deleter_ = nullptr;
  BaseAllocator* allocator_ = nullptr;
  size_t* shared_count_ = nullptr;
  size_t* weak_count_ = nullptr;
  //

  void destroy() {
    ptr_ = nullptr;
    deleter_ = nullptr;
    allocator_ = nullptr;
    shared_count_ = nullptr;
    weak_count_ = nullptr;
  }

  template<class U>
  void assigne(U&& other) {
    ptr_ = other.ptr_;
    deleter_ = other.deleter_;
    allocator_ = other.allocator_;
    shared_count_ = other.shared_count_;
    weak_count_ = other.weak_count_;
  }
 public:
  WeakPtr() = default;
  ~WeakPtr() {
    if (weak_count_ == nullptr) return;
    --*weak_count_;
    if (*weak_count_ + *shared_count_ == 0) {
      if (allocator_->TAllocated()) {
        allocator_->dealloc(ptr_);
      } else {
        allocator_->dealloc(deleter_);
      }
    }
  }

  size_t use_count() const {
    return *shared_count_;
  }

  bool expired() const {
    return !(*shared_count_);
  }

  SharedPtr<T> lock() const {
    if (expired())
      return SharedPtr<T>();
    return SharedPtr<T>(*this);
  }

  template<class U>
  WeakPtr(const SharedPtr<U>& other) {
    assigne(other);
    if (weak_count_)
      ++*weak_count_;
  }

  WeakPtr(const WeakPtr& other) {
    assigne(other);
    if (weak_count_)
      ++*weak_count_;
  }

  template<class U>
  WeakPtr(const WeakPtr<U>& other) {
    assigne(other);
    if (weak_count_)
      ++*weak_count_;
  }

  WeakPtr(WeakPtr&& other) {
    assigne(std::move(other));
    other.destroy();
  }

  template<class U>
  WeakPtr(WeakPtr<U>&& other) {
    assigne(std::move(other));
    other.destroy();
  }

  void swap(WeakPtr& other) {
    std::swap(ptr_, other.ptr_);
    std::swap(deleter_, other.deleter_);
    std::swap(allocator_, other.allocator_);
    std::swap(shared_count_, other.shared_count_);
    std::swap(weak_count_, other.weak_count_);
  }

  WeakPtr& operator=(const WeakPtr& other) {
    WeakPtr copy = other;
    swap(copy);
    return *this;
  }

  template<class U>
  WeakPtr& operator=(const WeakPtr<U>& other) {
    WeakPtr copy = other;
    swap(copy);
    return *this;
  }

  WeakPtr& operator=(WeakPtr&& other) {
    WeakPtr copy = std::move(other);
    swap(copy);
    return *this;
  }

  template<class U>
  WeakPtr& operator=(WeakPtr<U>&& other) {
    WeakPtr copy = std::move(other);
    swap(copy);
    return *this;
  }
};

