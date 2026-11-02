/* Progressive enhancement: all navigation works without JavaScript. */
const filter = document.getElementById('module-filter');
const modules = [...document.querySelectorAll('.module-list li')];
if (window.matchMedia('(max-width: 760px)').matches) {
  document.querySelector('.module-browser').open = false;
}
if (filter) {
  document.querySelector('.filter-controls').hidden = false;
  filter.addEventListener('input', () => {
    const query = filter.value.trim().toLowerCase();
    let visible = 0;
    for (const item of modules) {
      item.hidden = !item.textContent.toLowerCase().includes(query);
      if (!item.hidden) visible += 1;
    }
    document.getElementById('filter-status').textContent = `${visible} modules found`;
  });
}
