window.addEventListener('beforeunload', (event) => {
    fetch(window.location.origin + '/shutdown', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        keepalive: true,
        body: JSON.stringify({reason: 'User closed the window'})
    }).then(response => {
        console.log('Shutdown request sent.');
    }).catch(error => {
        console.error('Failed to send shutdown request:', error);
    });
});