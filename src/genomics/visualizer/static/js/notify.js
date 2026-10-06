// Desktop notifications when a background job finishes.
//
// Training and prediction jobs run for minutes to hours, so the tab is usually not the one in front
// when they end. A toast only helps if you are looking; this posts a system notification instead,
// but only when the tab is hidden (when it is visible the toast has already said it) and only after
// the user has asked for it — the browser requires a click to grant permission anyway.
const KEY = 'gv.notifyOnFinish';

export function notificationsSupported() {
  return typeof window !== 'undefined' && 'Notification' in window;
}

export function notificationsEnabled() {
  if (!notificationsSupported() || Notification.permission !== 'granted') return false;
  try { return localStorage.getItem(KEY) === '1'; } catch (e) { return false; }
}

/** Ask the browser for permission (must be called from a click); resolves to whether it is on. */
export async function setNotifications(on) {
  if (!notificationsSupported()) return false;
  if (!on) {
    try { localStorage.setItem(KEY, '0'); } catch (e) { /* private mode */ }
    return false;
  }
  let permission = Notification.permission;
  if (permission === 'default') {
    try { permission = await Notification.requestPermission(); } catch (e) { permission = 'denied'; }
  }
  const granted = permission === 'granted';
  try { localStorage.setItem(KEY, granted ? '1' : '0'); } catch (e) { /* private mode */ }
  return granted;
}

/** Blocked at the browser level: the page cannot re-ask, the user has to change it in site settings. */
export function notificationsBlocked() {
  return notificationsSupported() && Notification.permission === 'denied';
}

export function notifyFinished(title, body, { onClick } = {}) {
  if (!notificationsEnabled() || document.visibilityState === 'visible') return;
  try {
    const n = new Notification(title, { body, tag: `genomics-${title}` });
    n.onclick = () => { window.focus(); n.close(); if (onClick) onClick(); };
  } catch (e) { /* some browsers only allow notifications from a service worker */ }
}
